#!/usr/bin/env bash
## Deletes old tags for a Docker Hub repository, keeping "latest" plus the N
## most recently updated tags. Requires DOCKERHUB_USERNAME and DOCKERHUB_TOKEN
## in the environment (the token needs read/write/delete permissions).
set -euo pipefail

IMAGE="$1"        # e.g. myuser/myrepo-hail
KEEP="${2:-5}"    # number of most-recent non-"latest" tags to keep
PROTECT="${3:-}"  # tag this run just pushed (e.g. the short SHA) -- never delete it

TOKEN=$(curl -s -H "Content-Type: application/json" \
  -X POST -d "{\"username\": \"${DOCKERHUB_USERNAME}\", \"password\": \"${DOCKERHUB_TOKEN}\"}" \
  https://hub.docker.com/v2/users/login/ | jq -r .token)

if [ -z "$TOKEN" ] || [ "$TOKEN" == "null" ]; then
  echo "Failed to authenticate with the Docker Hub API" >&2
  exit 1
fi

# Docker Hub paginates tag listings, and only the newest tags (which we keep
# anyway) are needed to decide what to delete -- but we must walk every page
# to find and delete the old ones sitting past page 1.
TAGS=()
PAGE_URL="https://hub.docker.com/v2/repositories/${IMAGE}/tags/?page_size=100&ordering=-last_updated"
while [ -n "$PAGE_URL" ] && [ "$PAGE_URL" != "null" ]; do
  RESPONSE=$(curl -s -H "Authorization: JWT ${TOKEN}" "$PAGE_URL")
  while IFS= read -r name; do
    TAGS+=("$name")
  done < <(echo "$RESPONSE" | jq -r '.results[].name')
  PAGE_URL=$(echo "$RESPONSE" | jq -r '.next')
done

echo "Found ${#TAGS[@]} tags for ${IMAGE} (newest first)."

KEPT=0
for TAG in "${TAGS[@]}"; do
  if [ "$TAG" == "latest" ]; then
    continue
  fi
  # Docker Hub's ordering=-last_updated listing isn't guaranteed to reflect
  # a tag pushed moments ago (this script runs right after the push, in the
  # same job) -- it can sort as if older than every pre-existing tag, which
  # would otherwise get it treated as "oldest" and deleted. Protect this
  # run's own tag explicitly rather than trusting the API's freshness.
  if [ -n "$PROTECT" ] && [ "$TAG" == "$PROTECT" ]; then
    continue
  fi
  if [ "$KEPT" -lt "$KEEP" ]; then
    KEPT=$((KEPT + 1))
    continue
  fi
  echo "Deleting old tag: ${IMAGE}:${TAG}"
  curl -s -o /dev/null -w "  -> HTTP %{http_code}\n" -X DELETE \
    -H "Authorization: JWT ${TOKEN}" \
    "https://hub.docker.com/v2/repositories/${IMAGE}/tags/${TAG}/"
done

#!/bin/bash
set -eu

# Deploy documentation to GitHub Pages
# Usage: ./deploy-ghpages.sh [docs_directory]

DOCS_DIR="${1:-docs/html}"

if [ ! -d "$DOCS_DIR" ]; then
    echo "Error: Documentation directory '$DOCS_DIR' does not exist"
    exit 1
fi

echo "Deploying documentation from $DOCS_DIR to gh-pages branch..."

# Configure git
git config --global user.email "ci@circleci.com"
git config --global user.name "CircleCI"

# Create temporary directory for gh-pages
TMP_DIR=$(mktemp -d)
trap "rm -rf $TMP_DIR" EXIT

# Initialize a fresh repo with gh-pages as default branch
cd "$TMP_DIR"
git init -b gh-pages

# Copy documentation files
cp -r "$CIRCLE_WORKING_DIRECTORY/$DOCS_DIR/"* .

# Create .nojekyll to allow underscores in filenames
touch .nojekyll

# Commit all content
git add -A
git commit -m "Documentation update from ${CIRCLE_SHA1:0:7} on ${CIRCLE_BRANCH}"

# Set up remote with authentication token
REPO_URL="https://${GITHUB_TOKEN}@github.com/${CIRCLE_PROJECT_USERNAME}/${CIRCLE_PROJECT_REPONAME}.git"

# Force push to gh-pages branch (silencing output to avoid token exposure)
git push --force --quiet "$REPO_URL" gh-pages 2>&1 | grep -v "^remote:" || true

echo "Documentation successfully deployed to GitHub Pages!"

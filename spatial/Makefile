.PHONY: sync-core test deploy

CORE_REPO  ?= $(HOME)/Github/spatial-viewer-core
CORE_DIR   := output/viewer/core
DEPLOY_DIR := output/deploy

sync-core:
	@if [ ! -d "$(CORE_REPO)" ]; then \
	  echo "ERROR: $(CORE_REPO) not found. Clone it first:"; \
	  echo "  git clone https://github.com/stripathy/spatial-viewer-core $(CORE_REPO)"; \
	  exit 1; \
	fi
	@cd $(CORE_REPO) && git pull --ff-only
	mkdir -p $(CORE_DIR)
	rsync -av --delete \
	  --exclude='.git' --exclude='tests' --exclude='examples' \
	  --exclude='node_modules' --exclude='package*.json' --exclude='*.config.js' \
	  $(CORE_REPO)/src/ $(CORE_DIR)/
	@cd $(CORE_REPO) && git rev-parse --short HEAD > $(CURDIR)/$(CORE_DIR)/VERSION
	@echo ""
	@echo "Core synced. Pinned to $$(cat $(CORE_DIR)/VERSION)."
	@echo "Now: git diff $(CORE_DIR) and commit if the diff looks right."

test:
	npx playwright test

# SCZ deploys from output/deploy/ (which has its own copy of viewer files +
# per-sample JSON). Mirror viewer.js, index.html, style.css, scz-data.js,
# AND core/ before pushing to Netlify.
deploy:
	cp output/viewer/index.html output/viewer/style.css output/viewer/viewer.js output/viewer/scz-data.js $(DEPLOY_DIR)/
	mkdir -p $(DEPLOY_DIR)/core
	rsync -av --delete $(CORE_DIR)/ $(DEPLOY_DIR)/core/
	cd $(DEPLOY_DIR) && netlify deploy --prod --dir=.

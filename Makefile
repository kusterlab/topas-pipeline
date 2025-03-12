TOPAS_PIPELINE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(TOPAS_PIPELINE_DIR)MakefileShared

save_git_hash:
	git describe --dirty --always > hash.file

simsi: rm_err_file save_git_hash
	$(DOCKER_CMD) $(IMAGE) python3 -u -m bin.simsi -c $(LOCAL_DIR)/$(CONFIG_FILE) || (echo "1" > $(DATA)/err.out; exit 1)

full_pipeline: simsi save_git_hash
	$(DOCKER_CMD) $(IMAGE) python3 -u -m bin -c $(LOCAL_DIR)/$(CONFIG_FILE) || (echo "2" > $(DATA)/err.out; exit 2)

picked_group_fdr: simsi save_git_hash
	$(DOCKER_CMD) $(IMAGE) python3 -u -m bin.picked_group -c $(LOCAL_DIR)/$(CONFIG_FILE) || (echo "1" > $(DATA)/err.out; exit 1)

# runs sarcoma cohort in docker
docker_all: pull full_pipeline
docker_pgf: picked_group_fdr

# runs mintest
docker_mintest: CONFIG_FILE=config_patients_minimal_test.json
docker_mintest: full_pipeline

# runs glioma cohort in docker
docker_glioma: CONFIG_FILE=config_patients_glioma.json
docker_glioma: full_pipeline

# runs melanoma cohort in docker
docker_melanoma: CONFIG_FILE=config_melanoma.json
docker_melanoma: full_pipeline

# runs melanoma cohort in docker
docker_decryptm: CONFIG_FILE=config_patients_with_decryptm.json
docker_decryptm: full_pipeline


#plot scripts
heatmap: 
	python3 -u -m bin.heatmap_generator  -t plot_config.toml 


#genomics annotations
genomics_annotations:
	python3 -u -m bin.Oncostar_genomics_alterations && python3 -u -m bin.oncoKB_annotations 


# runs pipeline locally
all: DOCKER_CMD=
all: IMAGE=
all: LOCAL_DIR=.
all: full_pipeline

# runs glioma cohort locally
glioma: CONFIG_FILE=config_patients_glioma.json
glioma: DOCKER_CMD=
glioma: IMAGE=
glioma: LOCAL_DIR=.
glioma: full_pipeline

# runs melanoma cohort locally
melanoma: CONFIG_FILE=config_melanoma.json
melanoma: DOCKER_CMD=
melanoma: IMAGE=
melanoma: LOCAL_DIR=.
melanoma: full_pipeline

# runs minimal test for sarcoma cohort locally
mintest: CONFIG_FILE=config_patients_minimal_test.json
mintest: DOCKER_CMD=
mintest: IMAGE=
mintest: LOCAL_DIR=.
mintest: full_pipeline

# runs minimal test for models only cohort locally
models_only: CONFIG_FILE=config_patients_models_only_local.json
models_only: DOCKER_CMD=
models_only: IMAGE=
models_only: LOCAL_DIR=.
models_only: full_pipeline

# runs minimal test for glioma cohort locally
mintest_glioma: CONFIG_FILE=config_patients_glioma_minimal_test.json
mintest_glioma: DOCKER_CMD=
mintest_glioma: IMAGE=
mintest_glioma: LOCAL_DIR=.
mintest_glioma: full_pipeline

# runs minimal test for melanoma cohort locally
mintest_melanoma: CONFIG_FILE=config_melanoma_minimal_test.json
mintest_melanoma: DOCKER_CMD=
mintest_melanoma: IMAGE=
mintest_melanoma: LOCAL_DIR=.
mintest_melanoma: full_pipeline

test:
	python3 -m pytest --cov=bin --cov-report html --cov-report term tests/unit_tests


dependencies:
	git config --global credential.helper cache

build: dependencies save_git_hash
	docker build -f Dockerfile -t $(IMAGE) . || (exit 1)

registry:
	docker login gitlab.lrz.de:5005
	docker build -t gitlab.lrz.de:5005/proteomics/topas/topas-pipeline .
	docker push gitlab.lrz.de:5005/proteomics/topas/topas-pipeline

pull:
	docker login gitlab.lrz.de:5005
	docker pull gitlab.lrz.de:5005/proteomics/topas/topas-pipeline:master

jump:
	$(DOCKER_CMD) \
		$(IMAGE) bash

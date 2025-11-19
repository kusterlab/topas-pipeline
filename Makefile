TOPAS_PIPELINE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(TOPAS_PIPELINE_DIR)MakefileShared

save_git_hash:
	git describe --dirty --always > hash.file

full_pipeline: save_git_hash
	$(DOCKER_CMD) $(IMAGE) python3 -u -m topas_pipeline -c $(LOCAL_DIR)/$(CONFIG_FILE) || (echo "2" > $(DATA)/err.out; exit 2)

simsi: rm_err_file save_git_hash
	$(DOCKER_CMD) $(IMAGE) python3 -u -m topas_pipeline.simsi -c $(LOCAL_DIR)/$(CONFIG_FILE) || (echo "1" > $(DATA)/err.out; exit 1)

maxquant_qc: rm_err_file save_git_hash
	$(DOCKER_CMD) $(IMAGE) python3 -u -m topas_pipeline.search_qc.maxquant_qc -c $(LOCAL_DIR)/$(CONFIG_FILE) || (echo "1" > $(DATA)/err.out; exit 1)

picked_group_fdr: save_git_hash
	$(DOCKER_CMD) $(IMAGE) python3 -u -m topas_pipeline.picked_group -c $(LOCAL_DIR)/$(CONFIG_FILE) || (echo "1" > $(DATA)/err.out; exit 1)

report_creation: save_git_hash
	$(DOCKER_CMD) $(IMAGE) python3 -u -m topas_pipeline.report_creation -c $(LOCAL_DIR)/$(CONFIG_FILE) || (echo "1" > $(DATA)/err.out; exit 1)

# runs cohort in docker
docker_all: pull full_pipeline

# runs mintest
docker_mintest: CONFIG_FILE=config_patients_minimal_test.toml
docker_mintest: full_pipeline

# runs pipeline locally
all: DOCKER_CMD=
all: IMAGE=
all: LOCAL_DIR=.
all: full_pipeline

# runs minimal test locally
mintest: CONFIG_FILE=config_patients_minimal_test.toml
mintest: DOCKER_CMD=poetry run
mintest: IMAGE=
mintest: LOCAL_DIR=.
mintest: full_pipeline

#######################
## TESTING/PROFILING ##
#######################

test:
	python3 -m pytest --cov=topas_pipeline --cov-report html --cov-report term tests/unit_tests

mprof:
	mprof run --include-children --backend psutil_pss python3 -u -m topas_pipeline -c $(CONFIG_FILE)

mprof_plot:
	mprof plot -o mprofile_plot.png --backend agg


############
## DOCKER ##
############

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

bootstrap: DATA=/root/data
bootstrap:
	bash -c "cp /root/topas-pipeline/Makefile* $(LOCAL_DIR)"
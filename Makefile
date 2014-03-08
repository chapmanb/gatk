# gatk-framework -- Create executable script

all:
	rm -f target/*.jar
	rm -f public/gatk-package/target/*.jar
	mvn verify
	cat bin/gatk-framework.template target/GenomeAnalysisTK.jar > bin/gatk-framework
	chmod +x bin/gatk-framework

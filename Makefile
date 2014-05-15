# gatk-framework -- Create executable script

version=3.1-1

all:
	rm -f target/*.jar
	rm -f public/gatk-package/target/*.jar
	mvn verify
	rm -rf gatk-framework-${version}
	mkdir -p gatk-framework-${version}
	cp bin/gatk-framework gatk-framework-${version}
	chmod +x gatk-framework-${version}/gatk-framework
	cp target/GenomeAnalysisTK.jar gatk-framework-${version}/gatk-framework.jar
	tar -czvpf gatk-framework-${version}.tar.gz gatk-framework-${version}

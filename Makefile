# gatk-framework -- Create executable script

version_base=3.2
version=${version_base}-4

all:
	rm -f target/*.jar
	rm -f public/gatk-package-distribution/target/*.jar
	rm -f public/gatk-tools-public/target/*.jar
	rm -f public/gatk-utils/target/*.jar
	rm -f public/gatk-engine/target/*.jar
	mvn verify
	rm -rf gatk-framework-${version}
	mkdir -p gatk-framework-${version}
	cp bin/gatk-framework gatk-framework-${version}
	chmod +x gatk-framework-${version}/gatk-framework
	cp public/gatk-package-distribution/target/gatk-package-distribution-${version_base}*.jar gatk-framework-${version}/gatk-framework.jar
	tar -czvpf gatk-framework-${version}.tar.gz gatk-framework-${version}

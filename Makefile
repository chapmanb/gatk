# gatk-framework -- Create executable script

version_base=3.5
version=${version_base}-21

all:
	rm -f target/*.jar
	rm -f public/gatk-package-distribution/dependency-reduced-pom.xml
	rm -rf public/gatk-package-distribution/target
	rm -rf public/gatk-tools-public/target
	rm -rf public/gatk-utils/target
	rm -rf public/gatk-engine/target
	mvn clean
	mvn -e verify
	rm -rf gatk-framework-${version}
	mkdir -p gatk-framework-${version}
	cp bin/gatk-framework gatk-framework-${version}
	chmod +x gatk-framework-${version}/gatk-framework
	cp public/gatk-package-distribution/target/gatk-package-distribution-*.jar gatk-framework-${version}/gatk-framework.jar
	tar -czvpf gatk-framework-${version}.tar.gz gatk-framework-${version}

grep -v "^##" calls/popdel/polaris150.rnd.popdel.vcf | cut -f10- | grep -v "\./\." | sed -e 's;\([0-1]/[0-1]\)[[:graph:]]*;\1;g' -e 's;0/0;0;g' -e 's;0/1;1;g' -e 's;1/1;2;g' > calculations/pca/popdel.calls
#grep -v "^##" calls/delly/polaris150.delly.DEL.PASS.vcf | cut -f10- | grep -v "\./\." | sed -e 's;\([0-1]/[0-1]\)[[:graph:]]*;\1;g' -e 's;0/0;0;g' -e 's;0/1;1;g' -e 's;1/1;2;g' -e 's;-N1-DNA1-WGS1;;g'> calculations/pca/delly.calls
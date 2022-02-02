file(GLOB data *.root)
install( FILES ${data} DESTINATION ${product}/${version}/data)

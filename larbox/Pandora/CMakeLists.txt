install_fhicl()

file(GLOB xmls *.xml)
install( FILES ${xmls} DESTINATION ${product}/${version}/data)

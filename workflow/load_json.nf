import groovy.json.JsonSlurper

def load_json_file(json_file_path){
    def jsonSlurper = new JsonSlurper()

    // new File object from your JSON file

    def ConfigFile = new File(json_file_path)

    // load the text from the JSON

    String ConfigJSON = ConfigFile.text


    // create a dictionary object from the JSON text

    def myConfig = jsonSlurper.parseText(ConfigJSON)

    // to access values in the dict
    // log.info "Some val: ${myConfig.some_val}"

    return myConfig
}

// this function is based on the code found here:
// https://groups.google.com/g/nextflow/c/qzsORfO5CFU/m/pYh-tEWXAgAJ


//def map_json = new groovy.json.JsonSlurper().parse(new File('x.json'))
//map_json.each { key, value ->
//  println "$key : $value"
//}
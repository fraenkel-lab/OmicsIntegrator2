# Description of the files in visualization/src

1. create\_style\_code\_for\_visualization.py

    This Python code 
 - extracts min and max values of the 'ProteinChange' attribute from all nodes in a network JSON file (graph_json.json).
 - automatically generates the JavasScript code that programs the visual styles specific to the current network sample.
 - converts a network JSON file into a JavaScript code. 


2. style\_master.js

    This JavaScript code programs the visual styles that are common to all network samples.


3. visualize\_OmicsIntegrator\_results\_with\_CytoscapeJS\_lib.html

    This HTML code contains the JavaScript code that programs the network visualization of the OmicsIntegrator results 
    using the cytoscape.js library.
    It codes with the cose layout in cytoscape.js.
    The HTML code implements the "Save JSON" functionality.
    It also implements a legend for the node colors & border colors to be displayed with the network.


4. visualize\_user\_saved\_JSON\_file\_with\_CytoscapeJS\_lib.html

    This HTML code contains the JavaScript code that programs the network visualization of the user saved JSON file 
    using the cytoscape.js library.
    The HTML code implements the "Save JSON" functionality.
    It also implements a legend for the node colors & border colors to be displayed with the network.


5. run\_visualization\_for\_OmicsIntegrator\_results.sh

    This bash script runs the visualization for the OmicsIntegrator results.
    It creates a sample-specific directory visualize_results_<sample_name>/
    which contains the network graph.js and style.js specific to the current network sample, and a HTML file.
    

6. run\_visualization\_for\_OmicsIntegrator\_results.bat

    This script is a Windows version of run\_visualization\_for\_OmicsIntegrator\_results.sh.


7. run\_visualization\_for\_user\_saved\_JSON\_file.sh

    This bash script runs the visualization for the user saved JSON file.
    It creates a sample-specific directory visualize\_saved\_JSON\_<sample_name>/
    which contains the network graph.js and style.js specific to the current network sample, and a HTML file.


8. run\_visualization\_for\_user\_saved\_JSON\_file.bat

    This script is a Windows version of run\_visualization\_for\_user\_saved\_JSON\_file.sh.
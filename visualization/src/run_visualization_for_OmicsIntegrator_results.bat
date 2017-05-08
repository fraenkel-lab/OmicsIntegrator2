::
:: --------------------------------------------------------------
:: This script runs Visualization for the OmicsIntegrator results
:: --------------------------------------------------------------
::
SET FOREST_OUTPUT_DIR=%1
SET USER_SAMPLE=%2
SET USER_RESULT_DIR="../visualize_results_%USER_SAMPLE%"
SET CURRENT_DIR="%cd%"
::
:: Make a visualize_results_rundir/ after clearing out the existing one (if exists). 
:: This run dir is for creating the final graph.js and style.js files of a specific network, 
:: which will be used to visualize the network using the cytoscape.js library.
::
rmdir "../visualize_results_rundir" /s/q
mkdir "../visualize_results_rundir"
::
::
:: Concatenate the graph_head with the graph_json.json (which is the output from Forest).
::
cd %FOREST_OUTPUT_DIR%
copy "graph_json.json" %CURRENT_DIR%
cd %CURRENT_DIR%
copy /B "graph_head.txt" + "graph_json.json" "../visualize_results_rundir/graph.js"
::
::
:: Copy the style_master to the visualize_results_rundir/
::
copy "style_master.js" "../visualize_results_rundir/style.js"
::
::
:: Run Python code to extract min and max values from graph_json.json, to create the style.js code, 
:: and to process the graph.js.
::
py create_style_code_for_visualization.py
del "graph_json.json"
::
::
:: Now copy graph.js and style.js to a user-specific directory
:: and make the HTML file that contains JavaScript code using the cytoscape.js library to be user-specific in that directory,
:: so the HTML can be linked to a domain name and the user can access the weblink for their result
::
mkdir %USER_RESULT_DIR%
copy visualize_OmicsIntegrator_results_with_CytoscapeJS_lib.htm %USER_RESULT_DIR%
cd "../visualize_results_rundir"
copy graph.js %USER_RESULT_DIR%
copy style.js %USER_RESULT_DIR%
cd %USER_RESULT_DIR%
copy visualize_OmicsIntegrator_results_with_CytoscapeJS_lib.htm visualize_OmicsIntegrator_results_with_CytoscapeJS_lib_%USER_SAMPLE%.htm
del visualize_OmicsIntegrator_results_with_CytoscapeJS_lib.htm
cd %CURRENT_DIR%
::
::

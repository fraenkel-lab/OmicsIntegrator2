import json
import sys

def extract_min_max(json_file):
        """
        Extracts min and max values from a network JSON file

        Args:
                json_file: network JSON filename

        Returns:
                a tuple containing min and max values
        """

        json_data = open(json_file)
        network_data = json.load(json_data)
	
        node_list = network_data["elements"]["nodes"]
	
        min_value = 0
        max_value = 0
        for node in node_list:
            if "ProteinChange" in node["data"]:
                current_value = node["data"]["ProteinChange"]
            else:
                continue
            if current_value > max_value:
                max_value = current_value
            if current_value < min_value:
                min_value = current_value

        json_data.close()
	
        return [min_value, max_value]


def process_graph_json(graph_filename):
        """
        Converts a graph JSON file into a JavaScript code

        Args:
                graph_filename: graph filename

        Returns:
                None
        """
        
        graph_file = open(graph_filename, "a")
        graph_file.write(";")
        graph_file.close()
                        
        return None
        

def create_style_code_cytoscape(min_max, style_filename):
        """
        Creates the style JavaScript code to be used for the visualization
        using the cytoscape.js library.

        Args:
                min_max: a tuple containing min and max values
                style_filename: style filename

        Returns:
                None
        """

        min_value = min_max[0]
        max_value = min_max[1]
        min_string = str(min_value)
        max_string = str(max_value)
       
        style_file = open(style_filename, "a")

        style_file.write("\n")

        # Proteomic
        if max_value > 0:
            style_file.write("  }, {\n")
            style_file.write("    \"selector\" : \"node[ProteinChange > " + max_string + "]\",\n")
            style_file.write("    \"css\" : {\n")
            style_file.write("      \"background-blacken\" : 0.5\n")
            style_file.write("    }\n")
            style_file.write("  }, {\n")
            style_file.write("    \"selector\" : \"node[ProteinChange = " + max_string + "]\",\n")
            style_file.write("    \"css\" : {\n")
            style_file.write("      \"background-blacken\" : 0.5\n")
            style_file.write("    }\n")
            style_file.write("  }, {\n")
            style_file.write("    \"selector\" : \"node[ProteinChange > 0][ProteinChange < " + max_string + "]\",\n")
            style_file.write("    \"css\" : {\n")
            style_file.write("      \"background-blacken\" : \"mapData(ProteinChange, 0, " + max_string + ", -0.9, 0.5)\"\n")
            style_file.write("    }\n")
        if min_value < 0:
            style_file.write("  }, {\n")
            style_file.write("    \"selector\" : \"node[ProteinChange > " + min_string + "][ProteinChange < 0]\",\n")
            style_file.write("    \"css\" : {\n")
            style_file.write("      \"background-blacken\" : \"mapData(ProteinChange, " + min_string + ", 0, 0.5, -0.9)\"\n")
            style_file.write("    }\n")
            style_file.write("  }, {\n")
            style_file.write("    \"selector\" : \"node[ProteinChange = " + min_string + "]\",\n")
            style_file.write("    \"css\" : {\n")
            style_file.write("      \"background-blacken\" : 0.5\n")
            style_file.write("    }\n")
            style_file.write("  }, {\n")
            style_file.write("    \"selector\" : \"node[ProteinChange < " + min_string + "]\",\n")
            style_file.write("    \"css\" : {\n")
            style_file.write("      \"background-blacken\" : 0.5\n")
            style_file.write("    }\n")

        # Ending part
        style_file.write("  }, {\n")
        style_file.write("    \"selector\" : \"node:selected\",\n")
        style_file.write("    \"css\" : {\n")
        style_file.write("      \"background-color\" : \"rgb(255,255,0)\"\n")
        style_file.write("    }\n")
        style_file.write("  } ]\n")
        style_file.write("} ];")

        
        style_file.close()        
                      
        return None




if __name__ == '__main__':

        outdir = sys.argv[1]
        user_dir = sys.argv[2]

        process_graph_json("%s/graph.js"%outdir)

        min_max_values = extract_min_max("%s/graph_json.json"%user_dir)

        create_style_code_cytoscape(min_max_values, "%s/style.js"%outdir)


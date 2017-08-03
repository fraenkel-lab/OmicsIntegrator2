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
	
        min_pc = 0
        max_pc = 0
        for node in node_list:
            if "ProteinChange" in node["data"]:
                current_value = node["data"]["ProteinChange"]
            else:
                continue
            if current_value > max_pc:
                max_pc = current_value
            if current_value < min_pc:
                min_pc = current_value

        edge_list = network_data["elements"]["edges"]

        min_edge = 0
        max_edge = 0
        for edge in edge_list:
            current_value = edge["data"]["cost"]
            if current_value > max_edge:
                max_edge = current_value
            if current_value < min_edge:
                min_edge = current_value


        json_data.close()
	
        return [min_pc, max_pc, min_edge, max_edge]


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
        
        min_pc, max_pc, min_edge, max_edge = [str(item) for item in min_max]
       
        style_file = open(style_filename, "a")

        style_file.write("\n")

        # ProteinChange
        if max_pc > 0:
            style_file.write("  }, {\n"
                             "    \"selector\" : \"node[ProteinChange > " + max_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : 0.5\n"
                             "    }\n"
                             "  }, {\n"
                             "    \"selector\" : \"node[ProteinChange = " + max_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : 0.5\n"
                             "    }\n"
                             "  }, {\n"
                             "    \"selector\" : \"node[ProteinChange > 0][ProteinChange < " + max_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : \"mapData(ProteinChange, 0, " + max_pc + ", -0.9, 0.5)\"\n"
                             "    }\n")
        if min_pc < 0:
            style_file.write("  }, {\n"
                             "    \"selector\" : \"node[ProteinChange > " + min_pc + "][ProteinChange < 0]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : \"mapData(ProteinChange, " + min_pc + ", 0, 0.5, -0.9)\"\n"
                             "    }\n"
                             "  }, {\n"
                             "    \"selector\" : \"node[ProteinChange = " + min_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : 0.5\n"
                             "    }\n"
                             "  }, {\n"
                             "    \"selector\" : \"node[ProteinChange < " + min_pc + "]\",\n"
                             "    \"css\" : {\n"
                             "      \"background-blacken\" : 0.5\n"
                             "    }\n")
        # Edge Width
        style_file.write("  }, {\n"
                             "    \"selector\" : \"edge\",\n"
                             "  \"css\" : {\n"
                             "    \"width\" : \"mapData(cost," + max_edge + "," + min_edge + ",1,5)\"\n"
                             "  }\n")

        # Ending part
        style_file.write("  }, {\n"
                         "    \"selector\" : \"node:selected\",\n"
                         "    \"css\" : {\n"
                         "      \"background-color\" : \"rgb(255,255,0)\"\n"
                         "    }\n"
                         "  } ]\n"
                         "} ];")

        
        style_file.close()        
                      
        return None




if __name__ == '__main__':

        outdir = sys.argv[1]
        user_dir = sys.argv[2]

        process_graph_json("%s/graph.js"%outdir)

        min_max_values = extract_min_max("%s/graph_json.json"%user_dir)

        create_style_code_cytoscape(min_max_values, "%s/style.js"%outdir)


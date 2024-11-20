import pandas as pd
import networkx as nx
from .logger import *
import warnings

class ProteinInteractionNetwork:
    """
    Class to generate network for the predicted interactions.
    """

    def __init__(self, log = None):
        """
        Initialize ProteinInteractionNetwork class.

        :param log: Log file for logging messages.
        :type log: str
        """

        self.log = log
        
        return


    # Initiate graph
    def initiate_graph(self, interactions):
        """
        Create a NetworkX graph from a pandas DataFrame of interactions between hosts and pathogens.

        :param interactions: DataFrame containing predicted interactions between host and pathogen proteins.
        :type interactions: pandas.DataFrame

        :return: A graph object representing the interactions between hosts and pathogens.
        :rtype: networkx.Graph

        """

        # If no interactions
        if interactions.empty:
            warning_message = "No interactions found in this model - returning an empty graph"
            print(f"    > {warning_message}")
            self.log.warning(warning_message)
            return nx.Graph()
        else:
            try:
                G = nx.Graph()
                G.add_nodes_from(interactions['Host'], color="skyblue")
                G.add_nodes_from(interactions['Pathogen'], color="orange")
                G = nx.from_pandas_edgelist(interactions, 'Host', 'Pathogen', create_using=G)
                print("    > Graph object created")
                self.log.info("    > Graph object created")

                return G

            except Exception as e:
                print(f"HPIpy: Error occurred while creating the graph: {e}")
                self.log.error(f"Error occurred while creating the graph: {e}")
                raise
  

    # Protein hubs
    def calculate_protein_hubs(self, graph_object, networkOutDir):
        """
        Calculate protein hubs based on the degree centrality of nodes in a graph object,
        and store the results in an Excel file.

        :param graph_object: The graph object representing interactions between proteins.
        :param networkOutDir: The directory path where the output Excel file will be saved.
        
        :type graph_object: networkx.Graph
        :type networkOutDir: str

        :return: A DataFrame containing proteins and their degree centrality counts.
        :rtype: pd.DataFrame
        """

        if graph_object.number_of_nodes() == 0 or graph_object.number_of_edges() == 0:
            warning_message = "The graph is empty - unable to calculate protein hubs"
            print(f"    > {warning_message}")
            self.log.warning(warning_message)
            return {}
        else:
            try:
                protein_hubs_list = []

                for protein, degree in graph_object.degree():
                    protein_hubs_list.append({'Protein': protein, 'Count': degree})

                protein_hubs = pd.DataFrame(protein_hubs_list)
                protein_hubs.to_excel(f"{networkOutDir}/Protein_hubs.xlsx", index=False)
                print("    > Hubs calculation completed")
                self.log.info("    > Hubs calculation completed")

                return protein_hubs
        
            except Exception as e:
                print(f"HPIpy: Error calculating protein hubs: {e}")
                self.log.error(f"  Error calculating protein hubs: {e}")
                raise

    
    # Betweenness centrality
    def calculate_betweenness_centrality(self, graph_object, networkOutDir):
        """
        Calculate betweenness centrality for nodes in a graph and save the results to an Excel file.

        :param graph_object: NetworkX graph object for which betweenness centrality needs to be calculated.
        :param networkOutDir: Directory path where the Excel file will be saved.

        :type graph_object: networkx.Graph
        :type networkOutDir: str

        :return: Dictionary containing nodes as keys and their betweenness centrality as values.
        :rtype: dict
        """

        if graph_object.number_of_nodes() == 0 or graph_object.number_of_edges() == 0:
            warning_message = "The graph is empty - unable to calculate betweenness centrality"
            print(f"    > {warning_message}")
            self.log.warning(warning_message)
            return {}
        else:
            try:
                betweenness_centrality = nx.betweenness_centrality(graph_object, normalized=True, endpoints=True)
                betweenness_centrality_df = pd.DataFrame(betweenness_centrality.items(), columns=['Node', 'Betweenness_Centrality'])
                betweenness_centrality_df.to_excel(f'{networkOutDir}/Betweenness_centrality.xlsx', index=False)
                print("    > Betweenness centrality calculated")
                self.log.info("    > Betweenness centrality calculated")

                return betweenness_centrality
        
            except Exception as e:
                print(f"HPIpy: Error calculating betweenness centrality: {e}")
                self.log.error(f"Error calculating betweenness centrality: {e}")
                raise

    
    # Degree centrality
    def calculate_degree_centrality(self, graph_object, networkOutDir):
        """
        Calculate degree centrality for nodes in a graph and export the results to an Excel file.

        :param graph_object: NetworkX graph object for which degree centrality is to be calculated.
        :param networkOutDir: Directory path where the Excel file will be saved.

        :type graph_object: networkx.Graph
        :type networkOutDir: str

        :return: Dictionary of nodes with their degree centrality values.
        :rtype: dict
        """
        
        if graph_object.number_of_nodes() == 0 or graph_object.number_of_edges() == 0:
            warning_message = "The graph is empty - unable to calculate degree centrality"
            print(f"    > {warning_message}")
            self.log.warning(warning_message)
            return {}
        else:
            try:
                degree_centrality = nx.degree_centrality(graph_object)
                degree_centrality_df = pd.DataFrame(degree_centrality.items(), columns=['Node', 'Degree_Centrality'])
                degree_centrality_df.to_excel(f'{networkOutDir}/Degree_centrality.xlsx', index=False)
                print("    > Degree centrality calculated")
                self.log.info("    > Degree centrality calculated")
                
                return degree_centrality
        
            except Exception as e:
                print(f"HPIpy: Error calculating degree centrality: {e}")
                self.log.error(f"Error calculating degree centrality: {e}")
                raise

    
    # Closeness centrality
    def calculate_closeness_centrality(self, graph_object, networkOutDir):
        """
        Calculate closeness centrality for nodes in the graph and export the results to an Excel file.

        :param graph_object: The graph for which closeness centrality needs to be calculated.
        :param networkOutDir: Directory path where the Excel file will be saved.

        :type graph_object: networkx.Graph
        :type networkOutDir: str

        :return: Dictionary of closeness centrality values keyed by node.
        :rtype: dict

        :raises OSError: If the specified directory (networkOutDir) does not exist.
        """
        
        if graph_object.number_of_nodes() == 0 or graph_object.number_of_edges() == 0:
            warning_message = "The graph is empty - unable to calculate closeness centrality"
            print(f"    > {warning_message}")
            self.log.warning(warning_message)
            return {}
        else:
            try:
                closeness_centrality = nx.closeness_centrality(graph_object)
                closeness_centrality_df = pd.DataFrame(closeness_centrality.items(), columns=['Node', 'Closeness_Centrality'])
                closeness_centrality_df.to_excel(f'{networkOutDir}/Closeness_centrality.xlsx', index=False)
                print("    > Closeness centrality calculated")
                self.log.info("    > Closeness centrality calculated")

                return closeness_centrality
            
            except Exception as e:
                print(f"HPIpy: Error calculating closeness centrality: {e}")
                self.log.error(f"Error calculating closeness centrality: {e}")
                raise


    # Save network file
    def save_network(self, graph_object, networkOutDir):
        """
        Save a NetworkX graph object to a GML (Graph Modelling Language) file format.

        :param graph_object: NetworkX graph object to be saved.
        :param networkOutDir: Directory path where the GML file will be saved.
        
        :type graph_object: networkx.Graph
        :type networkOutDir: str

        :return: Path of the saved GML file.
        :rtype: str
        """

        if graph_object.number_of_nodes() == 0 or graph_object.number_of_edges() == 0:
            warning_message = "The graph is empty - no graph saved"
            print(f"    > {warning_message}")
            self.log.warning(warning_message)
            return {}
        else:
            try:
                save_graph = nx.write_gml(graph_object, f"{networkOutDir}/network.gml")
                print("    > Graph file saved successfully")
                self.log.info("    > Graph file saved successfully")

                return save_graph
            
            except Exception as e:
                print(f"HPIpy: Error saving network (GML) file: {e}")
                self.log.error(f"Error saving network (GML) file: {e}")
                raise

from .tables import *

class Interolog:
    """
    Class to perform interolog-based interaction prediction analysis between host and pathogen proteins.
    """

    def __init__(self):
        """
        Initialize the Interolog class.
        """

        self.connection = SQLconnect()


    # Filtering alignments
    def filter_blast(self, db, table, ident, evalue, cov, host = True):
        """
        Filters BLASTp results based on specific criteria.

        :param db: Name of the database.
        :param table: Table within the database to query.
        :param ident: Minimum identity threshold for filtering.
        :param evalue: Maximum E-value threshold for filtering.
        :param cov: Minimum coverage threshold for filtering.
        :param host: Whether to filter for host interactions (Default is True).

        :type db: str
        :type table: str
        :type ident: int
        :type evalue: float
        :type cov: int
        :type host: bool

        :return: A tuple containing the accession ID string and a DataFrame with filtered data.
        :rtype: tuple
        """
        
        eng = self.connection.create_connection(db)
        query = "SELECT * FROM {} WHERE pident >= {} AND evalue <= {} AND qcovs >= {} ".format(table, ident, evalue, cov)
        organisms = eng.execute(query).fetchall()
        
        accession_ids = []
        interaction_dict = {}
        
        for organism in organisms:
            unprocessed_id = organism[2]
            
            if unprocessed_id.startswith("sp") or unprocessed_id.startswith("tr") or unprocessed_id.startswith("ref"):
                processed_id = unprocessed_id.split("|")[1]
            elif unprocessed_id.startswith("dip"):
                processed_id = unprocessed_id.split("|")[0].split(":")[1]
            else:
                processed_id = unprocessed_id
            accession_ids.append(processed_id.split(".")[0])

            interaction_key = organism[1]
            interaction_value = processed_id.split(".")[0]

            if interaction_key in interaction_dict:
                interaction_dict[interaction_key].append(interaction_value)
            else:
                interaction_dict[interaction_key] = [interaction_value]

        p = []
        r = []
        vdf = pd.DataFrame()

        for key, values in interaction_dict.items():
            for value in values:
                p.append(key)
                r.append(value)

        if host:
            vdf = vdf.assign(Interactor_A=pd.Series(r))
            vdf = vdf.assign(Host=pd.Series(p))
            vdf.drop_duplicates(subset=None, keep='first', inplace=True)
        else:
            vdf = vdf.assign(Interactor_B=pd.Series(r, dtype = 'object'))
            vdf = vdf.assign(Pathogen=pd.Series(p, dtype = 'object'))
            vdf.drop_duplicates(subset=None, keep='first', inplace=True)

        accession_id_string = "(" + ",".join(f"'{accession_id}'" for accession_id in accession_ids) + ")"

        return accession_id_string, vdf


    # Searching the interactions
    def search_id(self, db, interactorA_filter, interactorB_filter, table):
        """
        Searches for interactinos between specified interactions in the database.

        :param db: Name of the database.
        :param interactorA_filter: Filter for Interactor A.
        :param interactorB_filter: Filter for Interactor B
        :param table: Table within the database to query.

        :type db: str
        :type interactorA_filter: str
        :type interactorB_filter: str
        :type table: str

        :return: DataFrame containing the interactions between specified interactors.
        :rtype: DataFrame
        """
        
        con = self.connection.create_connection(db)
        query = "SELECT * FROM {} WHERE Interactor_A IN {} AND Interactor_B IN {} OR Interactor_B IN {} AND Interactor_A IN {}".format(table, interactorA_filter, interactorB_filter, interactorA_filter, interactorB_filter)
        filtered = con.execute(query).fetchall()

        result = pd.DataFrame(filtered, columns = ['id', 'Interactor_A', 'Interactor_B', 'DetectionMethod', 'ConfidenceScore', 'PubMedID',
                                                   'EntrezGeneID_A', 'GO_IDs_A', 'GO_name_A', 'PDB_A', 'EntrezGeneID_B', 'GO_IDs_B', 'GO_name_B', 'PDB_B'])
        result = result.drop(columns=['id'])

        return result


    # To obtain the interactions
    def df_merge_interolog(self, result, pathogen_df, host_df):
        """
        Merges interaction results with host and pathogen DataFrames.

        :param result: DataFrame containing interaction results.
        :param pathogen_df: DataFrame containing pathogen data.
        :param host_df: DataFrame containing host data.
        
        :type result: DataFrame
        :type pathogen_df: DataFrame
        :type: host_df: DataFrame

        :return: Final DataFrame with merged interaction results.
        :rtype: DataFrame
        """
        
        merged_df = result.merge(pathogen_df, on='Interactor_B').merge(host_df, on='Interactor_A')
        final_df = merged_df[['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'DetectionMethod', 'ConfidenceScore', 'PubMedID',
                                                   'EntrezGeneID_A', 'GO_IDs_A', 'GO_name_A', 'PDB_A', 'EntrezGeneID_B', 'GO_IDs_B', 
                                                   'GO_name_B', 'PDB_B']].drop_duplicates()

        return final_df
    
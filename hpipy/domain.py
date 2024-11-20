from .tables import *

class Domain:
    """
    A class representing domain interactions and filtereing methods.
    """
    
    def __init__(self):
        """
        Initialize a Domain object with a SQL connection.
        """
        
        self.connection = SQLconnect()


    # Filtering domains
    def filter_hmm_domain(self, db, table, evalue, host=True):
        """
        Filter domain interactions based on hmmscan results.

        :param db: Name of the database.
        :param table: Table within the database to query.
        :param evalue: E-value threshold to filter hmmscan results.
        :param host: Whether to filter for host interactions (Defaults to True).
        
        :type db: str
        :type table: str
        :type evalue: float
        :type host: bool

        :return: A tuple containing a string of domain IDs and a DataFrame of interactions.
        :rtype: tuple
        """
        
        eng = self.connection.create_connection(db)
        query = "SELECT accessionT, query_name FROM {} WHERE E_valueB < {}".format(table, evalue)
        organism = eng.execute(query).fetchall()

        domain_accession = []
        domain_annotations = {}

        for record in organism:
            domain_accession.append(record[0])
            id = record[1]
            key = record[0]
            value = id
            if key in domain_annotations:
                domain_annotations[key].append(value)
            else:
                domain_annotations[key] = [value]

        interA = []
        interB = []
        domain_df = pd.DataFrame()
        for domain, annotations in domain_annotations.items():
            for annot in annotations:
                interA.append(domain.split(".")[0])
                interB.append(annot)

        if host:
            domain_df = domain_df.assign(Interactor_A=pd.Series(interA))
            domain_df = domain_df.assign(Host=pd.Series(interB))
            domain_df.drop_duplicates(subset=None, keep='first', inplace=True)
        else:
            domain_df = domain_df.assign(Interactor_B=pd.Series(interA))
            domain_df = domain_df.assign(Pathogen=pd.Series(interB))
            domain_df.drop_duplicates(subset=None, keep='first', inplace=True)
        
        domain_ids = [domain_id.split(".")[0] for domain_id in domain_accession]
        domain_id_st = "('" + "','".join(domain_ids) + "')"

        return domain_id_st, domain_df


    # Searching the interactions
    def search_domain_id(self, db, interactorA_filter, interactorB_filter, table):
        """
        Search for domain interactions in a database.

        :param db: Name of the database.
        :param interactorA_filter: Filter for Interactor A.
        :param interactorB_filter: Filter for Interactor B
        :param table: Table within the database to query.

        :type db: str
        :type interactorA_filter: str
        :type interactorB_filter: str
        :type table: str

        :return: A DataFrame containing the search results.
        :rtype: DataFrame
        """
        
        con = self.connection.create_connection(db)
        query = "SELECT * FROM {} WHERE Interactor_A IN {} AND Interactor_B IN {} OR Interactor_B IN {} AND Interactor_A IN {}".format(table, interactorA_filter, interactorB_filter, interactorA_filter, interactorB_filter)
        hmmscan_filter = con.execute(query).fetchall()

        result = pd.DataFrame(hmmscan_filter, columns = ['id', 'Interactor_A', 'Interactor_B', 'PfamName_A', 'InterproID_A', 'InterproName_A', 
                                                         'GO_IDs_A', 'GO_name_A', 'PDB_A', 'PfamName_B', 'InterproID_B', 'InterproName_B', 
                                                         'GO_IDs_B', 'GO_name_B', 'PDB_B'])
        result = result.drop(columns=['id'])

        return result
        

    # To obtain the interactions
    def df_merge_domain(self, result, pathogen_df, host_df):
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
        final_df = merged_df[['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'PfamName_A', 'InterproID_A', 'InterproName_A', 
                              'GO_IDs_A', 'GO_name_A', 'PDB_A', 'PfamName_B', 'InterproID_B', 'InterproName_B', 'GO_IDs_B', 
                              'GO_name_B', 'PDB_B']].drop_duplicates()

        return final_df
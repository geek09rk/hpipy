import sqlite3
from sqlite3 import Error
import pandas as pd

class SQLconnect:
    """
    Class for connecting to SQLlite databases.
    """

    def __init__(self):
        """
        Initialize the SQLconnect class.
        """
        
        return

    # Create connection to SQLite database
    def create_connection(self, db):
        """
        Create a connection to a SQLite database.

        :param db: The path to the SQLite database file.
        :type db: str

        :return: A connection object to the SQLite database.
        :rtype: sqlite3.Connection
        """
        
        conn = None
        try:
            conn = sqlite3.connect(db)
        except Error as e:
            print(e)
        
        return conn


class TABLES:
    """
    Class for creating tables in various databases.
    """

    def __init__(self, model = None, hostFile = None, pathogenFile = None):
        """
        Initialize the TABLES class.

        :param model: Model used for prediction.
        :param hostFile: The path to the host fasta file.
        :param pathogenFile: The path to the pathogen fasta file.

        :type model: str
        :type hostFile: str
        :type pathogenFile: str
        """
        
        self.model = model
        self.hostFile = hostFile
        self.pathogenFile = pathogenFile
        self.connection = SQLconnect()

        return
    

    # Tables for interaction databases
    def create_interaction_table(self, file = None, sep = None, table = None, db = None):
        """
        Create a table for interaction databases.

        :param file: The path to the file containing interaction data.
        :param sep: The delimiter used in the file.
        :param table: The name of the table to be created.
        :param db: The path to the SQLite database file.

        :type file: str
        :type sep: str
        :type table: str
        :type db: str

        :return: None
        """
        
        df = pd.read_csv(file, delimiter = sep)
        con = self.connection.create_connection(db)
        df.to_sql(table, con = con, if_exists = 'replace')
        con.close()

        return
    

    # BLAST results
    def create_table_blast(self, file = None, sep = None, table = None, db = None):
        """
        Create a table for BLASTp results in the database.

        :param file: The path to the file containing BLASTp results.
        :param sep: The delimiter used in the file.
        :param table: The name of the table to be created.
        :param db: The path to the SQLite database file.

        :type file: str
        :type sep: str
        :type table: str
        :type db: str

        :return: None
        """

        df = pd.read_csv(file, names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                                        'sstart', 'send', 'evalue', 'bitscore', 'qcovs'], delimiter = sep, low_memory = False)
        con = self.connection.create_connection(db)
        df.to_sql(table, con = con, if_exists = 'replace')
        con.close()

        return


    # hmmscan results
    def create_table_hmm(self, file = None, sep = None, table = None, db = None):
        """
        Create a table for hmmscan results in the database.

        :param file: The path to the file containing hmmscan results.
        :param sep: The delimiter used in the file.
        :param table: The name of the table to be created.
        :param db: The path to the SQLite database file.

        :type file: str
        :type sep: str
        :type table: str
        :type db: str

        :return: None
        """

        eng = self.connection.create_connection(db)
        df = pd.read_csv(file, names=['target_name', 'accessionT', 'query_name', 'accessionQ', 'E_valueF', 'scoreF',
                                    'biasF', 'E_valueB', 'scoreB', 'biasB', 'exp', 'reg', 'clu', 'ov', 'env', 'dom',
                                    'rep', 'inc', 'description_of_target'], delimiter = sep, comment = '#', index_col = False)
        df.to_sql(table, con = eng, if_exists = 'replace')
        eng.close()

        return
    

    # Human annotations table
    def table_human_annotations(self, file = None, sep = None, table = None, db = None):
        """
        Reads annotation file and store into a specified table in an SQLite database.

        :param file: The path to the CSV file to be read.
        :param sep: The delimiter used in the CSV file.
        :param table: The name of the table to be created in the SQLite database.
        :param db: The path to the SQLite database file.

        :type file: str
        :type sep: str
        :type table: str
        :type db: str

        :return: None
        """
        
        eng = self.connection.create_connection(db)
        df = pd.read_csv(file, delimiter = sep, index_col = False)
        df.to_sql(table, con = eng, if_exists = 'replace')
        eng.close()

        return
 

    # Query human protein annotations
    def query_human_annotations(self, proteinIDs, db, output_directory):
        """
        Queries an SQLite database for human annotations related to given protein IDs and saves the results to files.

        :param proteinIDs: A DataFrame containing the protein IDs to be queried. It should have a column 'Host' with the IDs.
        :param db: The path to the SQLite database file.
        :param output_directory: The directory where the query results will be saved as text files.
        
        :type proteinIDs: pandas.DataFrame
        :type db: str
        :type output_directory: str

        :return: A dictionary with table names as keys and the corresponding query results as pandas DataFrames.
        :rtype: dict
        """
        
        eng = self.connection.create_connection(db)
        
        tables = {
            'human_localization': ['Host', 'Localization'],
            'human_pathways': ['Host', 'Pathway', 'Description'],
            'human_drugs': ['Host', 'DrugID', 'Drug_common_name', 'GeneName', 'GenBankID'],
            'human_chembl': ['Host', 'ChEMBLID', 'ChEMBLName']
        }

        results = {}
        for table_name, columns in tables.items():
            ids = proteinIDs['Host'].apply(lambda x: x.split('|')[1] if '|' in x else x).tolist()

            table_results = []
            chunk_size = 1000
            
            for i in range(0, len(ids), chunk_size):
                chunk_ids = ids[i:i + chunk_size]
                
                query = f"""
                SELECT * FROM {table_name}
                WHERE Host IN ({','.join(['?']*len(chunk_ids))})
                """

                try:
                    df = pd.read_sql(query, con=eng, params=chunk_ids)
                    table_results.append(df)
                except Exception as e:
                    print(f"HPIpy: Error querying {table_name}: {e}")
            
            if table_results:
                combined_df = pd.concat(table_results, ignore_index=True)
                results[table_name] = combined_df
                combined_df.to_csv(f"{output_directory}/{table_name}.txt", sep="\t", index=False, columns=columns)

        eng.close()
        
        return results
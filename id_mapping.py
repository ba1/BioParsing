#!/home/bardya/anaconda3/envs/py3k/bin/python3

'''
Created on Oct 20, 2015

@author: bardya
'''
import os
import argparse
import sqlite3
import csv

def parse_args_subset1():
    parser = argparse.ArgumentParser(description='''Given a set of IDs of TYPE X in an input file in_ids.txt 
    find its corresponding id of TYPE Y and give it out else leave blank''')
    
    parser.add_argument('-mf', '--map-file', dest='mapfilepath', metavar='<file_path>', type=str, required=True,
                   help='path to file with accessions/ids map')
    
    parser.add_argument('-ms', '--map-file-separator', dest='mapfilesep', metavar='<str>', type=str, required=True, default='\t',
                   help='path to file with accessions/ids map')
    
    parser.add_argument('-i', '--input', dest='infilepath', metavar='<file_path>', type=argparse.FileType('rt'), required=True,
                   help='path to file with input accessions/ids')
    
    parser.add_argument('--version', action='version', version='0.1')
      
    parser.add_argument('-it', '--input-id-type', dest='intype', metavar='<int or str>', required=True,
                   help='''Select Column Number specifying the input type''')
    
    parser.add_argument('-ot', '--output-id-type', dest='outtype', metavar='<int or str>', required=True,
                   help='''Select Column Number specifying the output type''')
    
    return parser, parser.parse_args()


def parse_args_subset2(parent_parser, choices):
    parser = argparse.ArgumentParser(description='''Given a set of IDs of TYPE X in an input file in_ids.txt 
    find its corresponding id of TYPE Y and give it out else leave blank''', parents=[parent_parser], conflict_handler='resolve')

    parser.add_argument('-it', '--input-id-type', dest='intype', metavar='<int or str>', required=True, choices=choices,
                   help='''Select Column Number specifying the input type''')
    
    parser.add_argument('-ot', '--output-id-type', dest='outtype', metavar='<int or str>', required=True, choices=choices,
                   help='''Select Column Number specifying the output type''')

    
    return parser.parse_args()


def establishChoices(column_names):
    choices = []
    for i,nam in enumerate(column_names):
        if type(nam) is int:
            choices.extend([nam ,str(nam),i-len(column_names)])
        elif type(nam) is str:
            choices.extend([nam, nam.upper(), nam.lower()])
    
    return choices


def getColumnNames(map_path, sep='\t'):
    with open(map_path, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            else:
                return range(len(line.split(sep))+1)[1:]  #get a list starting from 1
        return range(len(line.split(sep))+1)[1:]  #get a list starting from 1

    
def map2sqldb(map_path, column_names, sep='\t'):
    """Determine the mean and 2std of the length distribution of a group
    """
    table_name = os.path.basename(map_path).rsplit('.', 1)[0]
    sqldb_name = table_name + '.sqlite3db'
    sqldb_path = os.path.join(os.path.dirname(map_path), sqldb_name)
    
    conn = sqlite3.connect(sqldb_path)  # @UndefinedVariable
    c = conn.cursor()

    # If table already exist, return the connector and the table_name
    SQL = '''
    SELECT count(*) FROM sqlite_master WHERE name == \"{}\"
    '''.format(table_name)
    c.execute(SQL)
    exists_flag = False
    if c.fetchone()[0] == 1:
        c.fetchall() #get rid of the remainder
        exists_flag=True
    
    if exists_flag:
        return c, table_name
    
        
    # Create table
    SQL = '''
    create table if not exists {0} ({1});
    '''.format(table_name, '\"' + '\" text,\"'.join([str(n).lower() for n in column_names]) + '\" text')
    
    c.execute(SQL)
    c.close()
    
    # Fill table
    SQL = '''
    insert into {0} values ({1})
    '''.format(table_name, ' ,'.join(['?']*len(column_names)))
    

    with open(map_path, 'r') as map_file:
        csv.field_size_limit(2147483647)
        csv_reader = csv.reader(map_file, delimiter=sep, quoting=csv.QUOTE_NONE)    
        
        with sqlite3.connect(sqldb_path) as conn:  # @UndefinedVariable
            c = conn.cursor()
            c.executemany(SQL, csv_reader)
    
    return c, table_name

       
def parseQuery(infilepath):
    query_ids = []
    with open(infilepath, 'r') as queryfile:
        for line in queryfile:
            query_ids.append(line.strip())
            
    return query_ids


def translateIDs(db_cursor, table_name, query_ids, intype, outtype):
    """Determine first column subject name for output file
    """
    target_ids = []
    
    for q in query_ids:
        mod_q = str(q)
        SQL = '''
        SELECT \"{}\" FROM {} WHERE \"{}\" LIKE \"{}\";
        '''.format(str(outtype), table_name, str(intype), mod_q)
        db_cursor.execute(SQL)
        target_ids.append(db_cursor.fetchall())
    
    return target_ids


def correctSeparator(sep):
    sepl = list(sep)
    if len(sepl) == 2 and sepl[0] == '\\':
        if sepl[1] == 'n':
            return '\n'
        if sepl[1] == 't':
            return '\t'
        return '\\' + sepl[1]
    else:
        return sep

if __name__ == '__main__':
    
    parent_parser, args1 = parse_args_subset1()
    
    sep = correctSeparator(args1.mapfilesep)
    column_names = getColumnNames(args1.mapfilepath, sep=sep)
    
    db_cursor, table_name = map2sqldb(args1.mapfilepath, column_names, sep=sep)
    query_ids = parseQuery(args1.infilepath.name)
    
    args2 = parse_args_subset2(parent_parser, establishChoices(column_names))
    target_ids = translateIDs(db_cursor, table_name, query_ids, args2.intype, args2.outtype)
    
    for t in target_ids:
        print(t)

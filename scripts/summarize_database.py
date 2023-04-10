import sys
sys.path.insert(1, '..')
sys.path.insert(1, '../src')
from src import database

data_handler = database.DatabaseHandler.create_cmod_handler()
disruptions_df = data_handler.query("select * from disruption_warning order by time",use_pandas=True)
grouped_df = disruptions_df.drop(columns=['commit_hash']).groupby(by=["shot"])
percent_nan = grouped_df.apply(lambda group: group.isnull().sum()*100/len(group))
percent_none = grouped_df.apply(lambda group: group.apply(lambda x: x == None, raw = True).sum()*100/len(group))
percent_nan.to_csv("../percent_nan.csv")
percent_none.to_csv("../percent_none.csv")
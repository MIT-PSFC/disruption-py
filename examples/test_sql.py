from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.settings.log_settings import LogSettings
from disruption_py.settings.output_type_request import OutputTypeRequest, ResultOutputTypeRequestParams, FinishOutputTypeRequestParams
from disruption_py.settings.shot_settings import ShotSettings
from disruption_py.utils.constants import TIME_CONST

cmod_handler = CModHandler()
cmod_database = cmod_handler.database

shots_to_check = cmod_database.query('select distinct shot from disruption_warning order by shot')['shot'].to_list()
# shots_to_check = [1140819003, 1140819002]

class LengthCheckOutputRequest(OutputTypeRequest):
    """
    Stream outputted data to a single csv file.
    
    Not recommended when retrieving a large number of shots. 
    """
    def __init__(self):
        self.disruptions = 0
        self.unequal_disruptions = 0
        self.non_disruptions = 0
        self.unequal_non_disruptions = 0

    def _output_shot(self, params : ResultOutputTypeRequestParams):
        sql_data = params.database.get_shots_data([params.shot_id], sql_table="disruption_warning")
        disrupted = params.result['time_until_disrupt'].notnull().all()
        
        if (
            len(sql_data) != len(params.result)
        ):
            min_diff = sql_data['time'].sort_values().diff().min()
        
        
        if disrupted:
            self.disruptions += 1
            new_format = (sql_data["dbkey"] > 1130000).all()
            
            # (new_format == True and len(sql_data) != len(params.result)) or \
            # (new_format == False and (len(sql_data) != len(params.result) + 1 and len(sql_data) != len(params.result)))
            if (
                # len(sql_data) != len(params.result) + 1 and 
                len(sql_data) != len(params.result)
            ):
                self.unequal_disruptions += 1
                print(f"Shot {params.shot_id} has {params.result.shape} shape but SQL has {len(sql_data)} rows for disruption. Min diff {min_diff}")
                
                
                if len(shots_to_check) < 5:
                    print("results:", params.result, "sql_data:", sql_data[['time', 'shot','dbkey', 'time_until_disrupt']])
        else: 
            self.non_disruptions += 1
            if len(sql_data) != len(params.result):
                self.unequal_non_disruptions += 1
                print(f"Shot {params.shot_id} has {params.result.shape} shape but SQL has {len(sql_data)} rows for no disruption. Min diff {min_diff}")
                if len(shots_to_check) < 5:
                    print("results:", params.result, "sql_data:", sql_data[['time', 'shot','dbkey', 'time_until_disrupt']])

    def get_results(self, params: FinishOutputTypeRequestParams):
        return self.disruptions, self.unequal_disruptions, self.non_disruptions, self.unequal_non_disruptions

shot_settings = ShotSettings(
    efit_tree_name="efit18",
    set_times_request="efit",
    run_columns=["time_until_disrupt"],
    run_tags=[],
    log_settings=LogSettings(
        log_to_console=False,
        log_file_path="examples/last_log.log",
    )
)
result_counts = cmod_handler.get_shots_data(
    shot_ids_request=shots_to_check,
    shot_settings=shot_settings,
    output_type_request=LengthCheckOutputRequest(),
    num_processes = 8,
)
print(result_counts)

# cmod_database.query('select * from disruption_warning where shot = 1120912003 ORDER BY time').to_csv('examples/data.csv')



"""
for shot_id in shots_to_check:
    sql_data = cmod_handler.database.get_shots_data([shot_id], sql_table="disruption_warning")
    sql_data_new = cmod_handler.database.get_shots_data([shot_id], sql_table="disruption_warning")
    if len(sql_data) != len(sql_data_new):
        print(f"Shot {shot_id} has {len(sql_data)} rows but SQL has {len(sql_data_new)} rows")

Shot 1140819002 has 752 rows but SQL has 188 rows
Shot 1140819003 has 760 rows but SQL has 190 rows
Shot 1140819004 has 768 rows but SQL has 192 rows
Shot 1150805012 has 0 rows but SQL has 84 rows
Shot 1150805013 has 0 rows but SQL has 89 rows
"""

"""
Shot 1140819001 has (97, 4) shape but SQL has 98 rows for no disruption
Shot 1140819003 has (95, 4) shape but SQL has 190 rows for no disruption
Shot 1140819002 has (94, 4) shape but SQL has 188 rows for no disruption
Shot 1140819004 has (96, 4) shape but SQL has 192 rows for no disruption
Shot 1150805013 has (88, 4) shape but SQL has 89 rows for no disruption
Shot 1150805014 has (85, 4) shape but SQL has 86 rows for no disruption
"""
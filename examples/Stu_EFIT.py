from disruption_py.machine.tokamak import Tokamak, resolve_tokamak_from_environment
from disruption_py.settings import RetrievalSettings
from disruption_py.workflow import get_shots_data


def main():
    """
    execute a simple workflow to fetch EFIT parameters.
    """

    tokamak = resolve_tokamak_from_environment()

    if tokamak is Tokamak.D3D:
        shotlist = [174960]
    else:
        raise ValueError(f"Unspecified or unsupported tokamak: {tokamak}.")

    print(f"Initialized for tokamak: {tokamak.value}")

    efit_params = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=RetrievalSettings(
        run_methods=["get_efit_parameters"],
        efit_nickname_setting="default"
        ),
        output_setting="dataset",
    )

    print(efit_params)

    geqdsk_params = get_shots_data(
        tokamak=tokamak,
        shotlist_setting=shotlist,
        retrieval_settings=RetrievalSettings(
            run_methods=["get_geqdsk_parameters"],
            efit_nickname_setting="default",
        ),
        output_setting="dataset",
    )

    print(geqdsk_params)

if __name__ == "__main__":
    main()
from Bio.PDB import PDBParser, PPBuilder
import plotly.express as px
import pandas as pd
from dash import Dash,dcc, html, Input, Output, no_update, ctx
import re
import argparse


def get_args():
    """
    Use to save the path of file as an argument in the terminal
    """
    parser = argparse.ArgumentParser(
        description="Take path of pdb file to use in the Ramachandran plot"
    )
    parser.add_argument("pdb_file_path", help="PDB file path, can be absolute or related to the current folder the script is in")

    return parser.parse_args()
def radian_to_degree_list(list_radian):
    """
    take for paramter a list of phi and psi angles in radian,
    and return a list of these angles in degree
    """
    list_degree= []
    for angles in list_radian:
        list_angles = []
        for angle in angles:
            if angle != None:
                list_angles.append(angle*57.2957795)
            elif angle ==  None:
                list_angles.append(angle)
        list_degree.append(list_angles)
    return list_degree


def pdb_dataframe(pdb_file):
    """
    take the path of a pdb file as a parameter,
    and return a dataframe with columns that contains selected informations from the pdb file
    """
    parser = PDBParser()
    id = re.search(".{4}[.]pdb$", pdb_file).group().rstrip(".pdb")
    protein = parser.get_structure(id , pdb_file)

    ppb = PPBuilder()
    model = protein[0]

    list_residue = []
    list_chain = []
    list_phi = []
    list_psi = []
    list_pp = []

    for chain in model:
        chain_id = str(chain).removeprefix("<").removesuffix(">")
        pps = ppb.build_peptides(chain)
        for pp in pps:
            phi_psi_list = radian_to_degree_list(pp.get_phi_psi_list())
            i = 0
            for residue in pp.get_sequence():
                list_phi.append(phi_psi_list[i][0])
                list_psi.append(phi_psi_list[i][1])
                list_residue.append(residue)
                list_chain.append(chain_id)
                list_pp.append(pp)
                i += 1

    df = pd.DataFrame([list_residue,list_phi,list_psi,list_chain])

    df_trans = df.transpose()
    df_trans.columns = ["residue", "phi", "psi","chain"]
    return df_trans


pdb_file = get_args().pdb_file_path
id = re.search(".{4}[.]pdb$",pdb_file).group().rstrip(".pdb")
df =  pdb_dataframe(pdb_file)
fig = px.scatter(df, x="phi", y="psi", hover_data=["residue", "chain"], color=None, range_x=[-180,180],
                 range_y=[-180,180],title=f"Ramachadran plot of {id}",width=600, height=600)


app = Dash(__name__)

app.layout = html.Div(children=[

    html.Button("No color", id="button_no_color"),
    html.Button("Color residues", id="button_residue"),
    html.Button("Color chains", id="button_chains"),

    dcc.Graph(id="Ramachadran plot", figure=fig)
])

@app.callback(
    Output("Ramachadran plot", "figure"),
    Input("button_no_color", "n_clicks"),
    Input("button_residue", "n_clicks"),
    Input("button_chains", "n_clicks")
)
def update(click1,click2,click3):

    if ctx.triggered_id == "button_no_color":

        fig = px.scatter(df, x="phi", y="psi", hover_data=["residue", "chain"], color=None, range_x=[-180, 180],
                         range_y=[-180, 180], title=f"Ramachadran plot of {id}",width=600, height=600)
        return fig

    if ctx.triggered_id == "button_residue":

        fig = px.scatter(df, x="phi", y="psi", hover_data=["residue", "chain"], color="residue", range_x=[-180, 180],
                         range_y=[-180, 180], title=f"Ramachadran plot of {id}",width=600, height=600)
        return fig

    if ctx.triggered_id == "button_chains":

        fig = px.scatter(df, x="phi", y="psi", hover_data=["residue", "chain"], color="chain", range_x=[-180, 180],
                         range_y=[-180, 180], title=f"Ramachadran plot of {id}",width=600, height=600)
        return fig

    return no_update


app.run()




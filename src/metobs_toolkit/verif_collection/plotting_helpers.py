
default_colorscheme = {
        #SCORES
        
        "RMSE": "#1f77b4",
        "MAE": "#ff7f0e",
        "MSE": "#2ca02c",
        "R2": "#d62728",
        "N": "#42027e",
        "modelbias": "#9467bd",
}





def get_color(column, colorscheme=default_colorscheme):
        if column in colorscheme:
                return colorscheme[column]
        return None


def create_linestyle_map(values):
        linestyle_list = ['-', '--', '-.', ':']
        linestyle_map = {val: linestyle_list[i % len(linestyle_list)] for i, val in enumerate(set(values))}
        return linestyle_map




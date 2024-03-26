#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

# Add title for each row in the legend
def add_row_titles( ax, num_rows, num_cols, row_title ):

    # original handles and labels
    handles, labels = ax.get_legend_handles_labels()

    # empty plot
    handles_title      = [ ax.plot( [], marker="", ls="" )[0] ]*(num_rows)

    # add the title
    handles_with_title = handles_title + handles
    labels_with_title  = row_title + labels

    # plot the legend
    legend = ax.legend(handles_with_title, labels_with_title, ncol=num_cols+1, loc='upper right', columnspacing=1.0)
    for vpack in legend._legend_handle_box.get_children()[:1]:
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)

# Add title for each column in the legend
def add_col_titles( ax, num_rows, num_cols, col_title ):

    # original handles and labels
    handles, labels = ax.get_legend_handles_labels()

    # empty plot
    handles_title      = [ ax.plot( [], marker="", ls="" )[0] ]*(num_cols)

    # add the title
    handles_with_title = []
    labels_with_title  = []
    for j in range(num_cols):
        handles_with_title += [handles_title[j]] + handles[j*num_rows:(j+1)*num_rows]
        labels_with_title  +=     [col_title[j]] +  labels[j*num_rows:(j+1)*num_rows]

    # plot the legend
    legend = ax.legend(handles_with_title, labels_with_title, ncol=num_cols, loc='upper right', columnspacing=1.0)
    for vpack in legend._legend_handle_box.get_children():
        for hpack in vpack.get_children()[:1]:
            hpack.get_children()[0].set_width(0)

def plot_example_figures_with_row_or_column_titles() -> None:

    # sampling points
    x = np.linspace(0.0, 2.0*np.pi, 100)

    # number of columns
    n_cols = 2

    # number of rows
    n_rows = 3

    # title of the rows and columns
    Row_Title = ["Row %d"%(row+1) for row in range(n_rows)]
    Col_Title = ["Col %d"%(col+1) for col in range(n_cols)]

    # create figure
    fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # plot the example curves, which are sine functions
    for j in range(n_cols):
        for i in range(n_rows):
            ax1.plot(x, np.sin(x-0.5*np.pi*j)-i, label='(%d, %d)'%(i+1, j+1))
            ax2.plot(x, np.sin(x-0.5*np.pi*j)-i, label='(%d, %d)'%(i+1, j+1))

    # add the title in the legend
    add_row_titles( ax1, n_rows, n_cols, Row_Title )
    add_col_titles( ax2, n_rows, n_cols, Col_Title )

    # save to file
    plt.tight_layout()
    fig.savefig( 'fig_legends_with_row_or_column_titles.png' )

def main() -> None:
    plot_example_figures_with_row_or_column_titles()

if __name__ == '__main__':
    main()

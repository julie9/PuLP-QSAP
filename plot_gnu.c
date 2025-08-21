#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void generate_gnuplot_file(int n, int p, int* assignment, int* timeOfTask) {
    /* Create a 2D array [p x n], initialized to 0. row -> processor, col -> task */
    int** processorData = (int**)malloc(p * sizeof(int*));
    for (int i = 0; i < p; i++) {
        processorData[i] = (int*)calloc(n, sizeof(int));
    }

    /* Fill processorData based on assignment/timeOfTask */
    for (int task = 0; task < n; task++) {
        int procIndex = assignment[task];
        if (procIndex >= 0 && procIndex < p) {
            processorData[procIndex][task] = timeOfTask[task];
        }
    }

    /* Write out the data file "data.dat" */
    FILE* dataFile = fopen("data.dat", "w");
    if (!dataFile) {
        fprintf(stderr, "Error: Could not open data.dat for writing.\n");
        return;
    }

    /* Each row:  label col1 col2 ... coln */
    for (int proc = 0; proc < p; proc++) {
        fprintf(dataFile, "p%d", proc);
        for (int task = 0; task < n; task++) {
            fprintf(dataFile, " %d", processorData[proc][task]);
        }
        fprintf(dataFile, "\n");
    }

    fclose(dataFile);

    /* Write out the gnuplot script "plot.gp". */
    FILE* gpFile = fopen("plot.gp", "w");
    if (!gpFile) {
        fprintf(stderr, "Error: Could not open plot.gp for writing.\n");
        return;
    }

    fprintf(gpFile,
        "set terminal pngcairo size 800,600\n"
        "set output 'plot.png'\n"
        "set title 'Stacked Bar Plot of Task Times per Processor'\n"
        "set style data histograms\n"
        "set style histogram rowstacked\n"
        "set style fill solid 1.0 border -1\n"
        "set auto x\n"
        "set grid y\n"
        "set boxwidth 0.75\n"
        "set xtics rotate by -45 font ',6'\n"  // Set font size to 8 for x-axis labels
        "unset key\n"  // Add this line to remove the legend
        "set palette defined (0 'red', 1 'green', 2 'blue', 3 'yellow', 4 'cyan', 5 'magenta')\n"  // Define a custom color palette
        "\n"
        "# Now we plot each task's column on top of each other:\n"
        "plot \\\n"
    );

    for (int task = 0; task < n; task++) {
        if (task == 0) {
            fprintf(gpFile, "  'data.dat' using %d:xtic(1) title 'Task%d' linecolor palette frac %f", task + 2, task, (float)task / n);
        } else {
            fprintf(gpFile, ", \\\n  '' using %d title 'Task%d' linecolor palette frac %f", task + 2, task, (float)task / n);
        }
    }

    fprintf(gpFile, "\n");

    fclose(gpFile);

    /* Free memory */
    for (int i = 0; i < p; i++) {
        free(processorData[i]);
    }
    free(processorData);

    printf("Data file 'data.dat' and gnuplot script 'plot.gp' have been created.\n");
    printf("Run the following to generate the plot:\n\n");
    printf("    gnuplot plot.gp\n\n");
    printf("This will produce 'plot.png' with the stacked bar chart.\n");
}

// int main(void) {
//     int n = 5;  /* number of tasks */
//     int p = 2;  /* number of processors */

//     /* Hardcoded example arrays */
//     int assignment[] = {0, 1, 0, 0, 1}; /* which processor each task goes to */
//     int timeOfTask[] = {1, 1, 1, 1, 1}; /* time for each task */

//     generate_gnuplot_file(n, p, assignment, timeOfTask);

//     return 0;
// }

// To compile and run:
// gcc -o plot_gnu plot_gnu.c
// ./plot_gnu
// gnuplot plot.gp
// This will generate a plot.png file with the stacked bar chart.
// The plot.png file will show a stacked bar chart with each task represented as a colored bar on the y-axis (processor) and the x-axis (time).
// The color of each bar will be based on the task number, with the color palette defined in the gnuplot script.
// The x-axis labels will be rotated by -45 degrees to make them more readable.
// The plot will be saved as a PNG file named plot.png.
// The gnuplot script plot.gp will contain the commands to generate the plot.
// The data.dat file will contain the data for the plot, with each row representing a processor and each column representing a task.
// The generate_gnuplot_file function takes the number of tasks, number of processors, an array of task assignments to processors, and an array of task times as input.
// The function creates a 2D array to store the data for the plot, fills the array based on the input data, and writes the data to a data.dat file.


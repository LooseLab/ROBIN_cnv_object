from cnv_from_bam import iterate_bam_file
from cnsmeth.subpages.base_analysis import BaseAnalysis
import natsort
from cnsmeth import theme, resources
import pandas as pd
import logging
import numpy as np
import os
import sys
import asyncio
from nicegui import ui, app
import click
from pathlib import Path
import pickle
import ruptures as rpt
from cbseg import determine_cbs_stat, determine_t_stat, determine_cbs, segment, validate
import statistics
from composition_stats import clr
import io
import tempfile

os.environ["CI"] = "1"


# Calculate average values
class Result:
    def __init__(self, cnv_dict):
        self.cnv = cnv_dict


def iterate_bam(
        bamfile, _threads=1, mapq_filter=60, copy_numbers=None, log_level=int(logging.ERROR)
):
    """
    Iterate over a BAM file and perform CNV analysis.

    Parameters:
        bamfile (str): Path to the BAM file.
        _threads (int, optional): Number of threads for analysis. Default 1.
        mapq_filter (int, optional): Mapping quality filter. Default 60.
        copy_numbers (list, optional): List of copy numbers. Default None.
        log_level (int, optional): Logging level. Default int(logging.ERROR).

    Returns:
        tuple: Tuple containing the result and copy numbers.
    """
    result = iterate_bam_file(
        bamfile,
        _threads=_threads,
        mapq_filter=mapq_filter,
        copy_numbers=copy_numbers,
        log_level=log_level,
    )
    return result, copy_numbers


def reduce_list(lst, max_length=1000):
    """
    Reduce the length of a list by half until it reaches the specified maximum length.

    Parameters:
        lst (list): Input list.
        max_length (int, optional): Maximum length of the list. Default 1000.

    Returns:
        list: Reduced list.
    """
    while len(lst) > max_length:
        lst = lst[::2]
    return lst


class CNV_Difference:
    def __init__(self, *args, **kwargs):
        self.cnv = {}


def moving_average(data, n=3):
    return np.convolve(data, np.ones(n) / n, mode="same")


class CNVAnalysis(BaseAnalysis):
    def __init__(self, *args, target_panel=None, **kwargs):
        self.penalty = None
        self.file_list = []
        self.cnv_dict = {}
        self.cnv_dict["bin_width"] = 0
        self.cnv_dict["variance"] = 0
        self.update_cnv_dict = {}
        self.result = None
        self.result2 = None
        self.result3 = CNV_Difference()
        self.element_counter = 1
        self.XYestimate = "Unknown"
        self.median_cnv_list = []
        self.changepoint_list_scaled = []
        self.changepoint_y = []
        self.contig_list = []
        self.changepoint_report = {}

        # Separate changepoint dictionaries
        self.changepoint_dict = {}
        self.reference_changepoint_dict = {}
        self.difference_changepoint_dict = {}
        self.changepoint_toggle = None
        self.step_line_toggle = None
        self.scatter_toggle = None

        with open(
                os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)),
                    "HG01280_control.pkl",
                ),
                "rb",
        ) as f:
            self.ref_cnv_dict = pickle.load(f)
        self.target_panel = target_panel
        if self.target_panel == "rCNS2":
            self.gene_bed_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "rCNS2_panel_name_uniq.bed",
            )
        elif self.target_panel == "AML":
            self.gene_bed_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "AML_panel_name_uniq.bed",
            )
        self.gene_bed = pd.read_table(
            os.path.join(self.gene_bed_file),
            names=["chrom", "start_pos", "end_pos", "gene"],
            header=None,
            sep="\s+",
        )
        super().__init__(*args, **kwargs)

    def generate_csv_report(self):
        # Create a DataFrame from the changepoint_dict
        report_data = []
        for chromosome, changepoints in self.changepoint_dict.items():
            for changepoint in changepoints:
                report_data.append({'Chromosome': chromosome, 'Changepoint': changepoint * self.cnv_dict["bin_width"]})

        df = pd.DataFrame(report_data)

        # Convert DataFrame to CSV
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)

        return csv_buffer.getvalue()

    def detect_changepoints(self, contig, cnv, penalty, total, changepoint_dict):
        changepoint_list = [total]
        changepoint_list_normal = []
        np_cnv = np.array(cnv)
        if contig == "chrM":
            return changepoint_list, changepoint_dict
        if len(np_cnv) > 3:
            algo = rpt.KernelCPD(kernel="linear", min_size=2).fit(np_cnv)
            changepoints = algo.predict(pen=int(penalty))
            for changepoint in changepoints:
                changepoint_list.append(changepoint + total)
                changepoint_list_normal.append(changepoint)
            changepoint_dict[contig] = changepoint_list_normal
        return changepoint_list, changepoint_dict

    def estimate_XY(self):
        # We remove zero points as they are likely centromeric.
        X = round(np.average([i for i in self.result3.cnv["chrX"] if i != 0]), 2)
        Y = round(np.average([i for i in self.result3.cnv["chrY"] if i != 0]), 2)
        # print (f"X={X},Y={Y}")
        if X >= 0.1 and Y <= 0.1:
            self.XYestimate = "XX"
        elif X <= 0.1 and Y >= -0.2:
            self.XYestimate = "XY"
        else:
            self.XYestimate = "Unknown"

    async def process_bam(self, bamfile, timestamp):
        self.file_list.append(bamfile)
        # cnv_dict = self.update_cnv_dict.copy()
        # self.result, self.update_cnv_dict = await run.cpu_bound(iterate_bam, bamfile, _threads=self.threads, mapq_filter=60, copy_numbers=cnv_dict)
        # print (f"Processing {bamfile}, {timestamp}")
        await self.do_cnv_work(bamfile)

    async def do_cnv_work(self, bamfile):
        # self.result, self.update_cnv_dict = background_tasks.create(run.cpu_bound(iterate_bam, bamfile, _threads=self.threads, mapq_filter=60, copy_numbers=self.update_cnv_dict))

        self.result = iterate_bam_file(
            bamfile,
            _threads=self.threads,
            mapq_filter=60,
            copy_numbers=self.update_cnv_dict,
            log_level=int(logging.ERROR),
        )

        self.cnv_dict["bin_width"] = self.result.bin_width
        self.cnv_dict["variance"] = self.result.variance
        self.result2 = iterate_bam_file(
            bam_file_path=None,
            _threads=self.threads,
            mapq_filter=60,
            copy_numbers=self.ref_cnv_dict,
            log_level=int(logging.ERROR),
            bin_width=self.cnv_dict["bin_width"],
        )

        for key in self.result.cnv.keys():
            if key != "chrM":
                # print(key, np.mean(self.result.cnv[key]))#[i for i in self.result.cnv[key] if i !=0]))
                moving_avg_data1 = moving_average(self.result.cnv[key])
                moving_avg_data2 = moving_average(self.result2.cnv[key])
                self.result3.cnv[key] = moving_avg_data1 - moving_avg_data2

        self.estimate_XY()

        if self.summary:
            with self.summary:
                self.summary.clear()
                with ui.row():
                    if self.XYestimate != "Unknown":
                        if self.XYestimate == "XY":
                            ui.icon("man").classes("text-4xl")
                        else:
                            ui.icon("woman").classes("text-4xl")
                        ui.label(f"Estimated Genetic Sex: {self.XYestimate}")
                    ui.label(f"Current Bin Width: {self.result.bin_width}")
                    ui.label(f"Current Variance: {round(self.result.variance, 3)}")
        np.save(os.path.join(self.output, "CNV.npy"), self.result.cnv)
        np.save(os.path.join(self.output, "CNV_dict.npy"), self.cnv_dict)

        # Only update the plot if the queue is empty?
        if self.bamqueue.empty() or self.bam_processed % 5 == 0:
            self._update_cnv_plot(
                plot_to_update=self.scatter_echart, result=self.result, title="CNV"
            )
            self._update_cnv_plot(
                plot_to_update=self.reference_scatter_echart,
                result=self.result2,
                title="Reference CNV",
            )
            self._update_cnv_plot(
                plot_to_update=self.difference_scatter_echart,
                result=self.result3,
                title="Difference CNV",
                min="dataMin",
            )
        # else:
        await asyncio.sleep(0.05)
        self.running = False

    def update_plots(self, gene_target=None, penalty=10):
        self.selected_contig = self.chrom_select.value  # Update the selected contig
        self._update_cnv_plot(
            plot_to_update=self.scatter_echart,
            result=self.result,
            gene_target=gene_target,
            penalty=penalty,
            title="CNV",
        )
        self._update_cnv_plot(
            plot_to_update=self.reference_scatter_echart,
            result=self.result2,
            gene_target=gene_target,
            penalty=penalty,
            title="Reference CNV",
        )
        self._update_cnv_plot(
            plot_to_update=self.difference_scatter_echart,
            result=self.result3,
            gene_target=gene_target,
            title="Difference CNV",
            penalty=penalty,
            min="dataMin",
        )

    def setup_ui(self):
        self.display_row = ui.row()
        if self.summary:
            with self.summary:
                ui.label("No CNV data available.")
        with self.display_row:
            # self.progrock.visible = False
            ui.label("Copy Number Variation").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
        with ui.row():
            self.pen_select = ui.number(
                on_change=lambda e: (self.update_plots(penalty=e.value)),
                label="Select Penalty",
                value=10,
                placeholder=10,
                min=0
            ).style("width: 100px")

            self.chrom_select = ui.select(
                options={"All": "All"},
                on_change=self.update_plots,
                label="Select Chromosome",
                value="All",
            ).style("width: 150px")
            self.gene_select = ui.select(
                options={"All": "All"},
                on_change=lambda e: (
                    self.update_plots()
                    if e.value == "All"
                    else self.update_plots(gene_target=e.value)
                ),
                label="Select Gene",
                value="All",
            ).style("width: 150px")
            ui.label().bind_text_from(
                self.cnv_dict, "bin_width", backward=lambda n: f"Bin Width: {n}"
            )
            ui.label().bind_text_from(
                self.cnv_dict, "variance", backward=lambda n: f"Variance: {round(n, 3)}"
            )
            # Toggle for step line
            self.step_line_toggle = ui.toggle(
                on_change=lambda e: self.update_plots(),
                options={"on": "Show Step Line", "off": "Hide Step Line"},
                value="on"
            ).style("width: 250px")
            # Toggle for rupture changepoints
            self.changepoint_toggle = ui.toggle(
                on_change=lambda e: self.update_plots(),
                options={"on": "Show Changepoints", "off": "Hide Cheangepoints"},
                value="off"
            ).style("width: 300px")

            # Toggle for scatter plot
            self.scatter_toggle = ui.toggle(
                on_change=lambda e: self.update_plots(),
                options={"on": "Show Scatter", "off": "Hide Scatter"},
                value="on"
            ).style("width: 220px")

            # Button for csv report download
            with ui.button(on_click=self.download_report):
                ui.image('download_sprite.png').classes('rounded-full w-20 h-20 ml-4')



        self.scatter_echart = self.generate_chart(title="CNV Scatter Plot")

        self.difference_scatter_echart = self.generate_chart(
            title="Difference Plot", initmin=-2, initmax=2
        )  # , type="log")

        with ui.expansion("See Reference DataSet", icon="loupe").classes("w-full"):
            self.reference_scatter_echart = self.generate_chart(
                title="Reference CNV Scatter Plot"
            )

    async def download_report(self):
        csv_data = self.generate_csv_report()

        # Use a temporary file to save the CSV data
        with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp_file:
            tmp_file.write(csv_data.encode('utf-8'))
            tmp_file_path = tmp_file.name

        # Trigger the download
        ui.download(tmp_file_path, 'changepoint_report.csv')


    def generate_chart(self, title=None, initmax=8, initmin=0, type="value"):
        return (
            ui.echart(
                {
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": f"{title}"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {
                        "type": f"{type}",
                        "max": "dataMax",
                        "splitLine": {"show": False},
                    },
                    "yAxis": {
                        "type": "value",
                        "logBase": 2,
                    },
                    "dataZoom": [
                        {"type": "slider", "filterMode": "none"},
                        {
                            "type": "slider",
                            "yAxisIndex": 0,
                            "filterMode": "none",
                            "startValue": initmin,
                            "endValue": initmax,
                        },
                        # {"type": "inside", "xAxisIndex": 0, "filterMode": "none"},
                        # {"type": "inside", "yAxisIndex": 0, "filterMode": "none"},
                    ],
                    "series": [
                        {
                            "type": "scatter",
                            "symbolSize": 1,
                            "data": [],
                        }
                    ],
                }
            )
            .style("height: 450px")
            .classes("border-double")
        )

    def _update_cnv_plot(
            self, plot_to_update=None, result=None, gene_target=None, title=None, min=0, penalty=10
    ):
        if penalty <= 0:
            ui.notify('The penalty must be positive!', type='error', color='red')
            return

        if result or self.result:
            cosmic_panel = pd.read_csv('COSMIC_census.csv')
            df = cosmic_panel["Genome Location"]
            # Split the 'genome_locations' column
            df1 = df.str.split(':', expand=True)
            df1.columns = ["chromosome", "position"]
            df2 = df1["position"].str.split('-', expand=True)
            df2.columns = ["start_pos", "end_pos"]
            # print(df2)
            frames = df1["chromosome"], df2
            # final_df = pd.join(frames)
            # print(final_df)

            total = 0
            valueslist = {"All": "All"}  # A dictionary to store chromosome data
            genevalueslist = {"All": "All"}  # A dictionary to store gene data
            self.chrom_filter = self.chrom_select.value  # chrom_filter = integer relating to the chromosome.

            # Select the appropriate changepoint dictionary
            if plot_to_update == self.scatter_echart:
                changepoint_dict = self.changepoint_dict
            elif plot_to_update == self.reference_scatter_echart:
                changepoint_dict = self.reference_changepoint_dict
            elif plot_to_update == self.difference_scatter_echart:
                changepoint_dict = self.difference_changepoint_dict

            # Min and max values for x-axis
            min = min
            max = "dataMax"
            # If gene_target is provided, retrieve start and end positions of the gene
            if gene_target:
                start_pos = self.gene_bed.iloc[int(gene_target)].start_pos
                end_pos = self.gene_bed.iloc[int(gene_target)].end_pos
                chrom = self.gene_bed.iloc[int(gene_target)].chrom
                for counter, contig in enumerate(
                        natsort.natsorted(result.cnv), start=1
                ):
                    valueslist[counter] = contig  # Key of valuelist is counter, contig is value
                    if contig == chrom:
                        break
                self.chrom_filter = counter  # Set chrom_filter based on counter
                min = start_pos - 10 * self.cnv_dict["bin_width"]  # Set min x value as 10 before start position
                max = end_pos + 10 * self.cnv_dict["bin_width"]  # Set max x value as 10 after end position
            if self.chrom_filter == "All":
                counter = 0

                plot_to_update.options["title"]["text"] = f"{title} - All Chromosomes"  # Update plot title text
                plot_to_update.options["series"] = []  # Clear series data
                cp_total = 0
                offset = 0
                # Modified for statement
                np_cnv_all = []
                median_chrom_cnv = []
                for contig, cnv in result.cnv.items():
                    # Reformat into numpy array for compatibility with ruptures.
                    np_cnv = np.array(cnv)
                    median_chrom_cnv.append(statistics.median(np_cnv))
                    np_cnv_all.append(np_cnv)
                combined_array = np.concatenate(np_cnv_all, axis=0)
                for contig, cnv in result.cnv.items():
                    np_cnv = np.array(cnv)

                    ## Log Scaling
                    # np_cnv = np_cnv + 0.001
                    # np_cnv = clr(np_cnv)
                    # print(clr_cnv)

                    # Skip mitochondria chrom if present
                    if contig == "chrM":
                        continue
                    counter += 1
                    valueslist[counter] = contig  # Add chromosome to valueslist
                    # Change point detection, provided chromosome consists of more than 3 data points.
                    if len(np_cnv) > 3:
                        changepoint_list, changepoint_dict = self.detect_changepoints(contig, cnv, penalty, total,
                                                                                      changepoint_dict)
                    offset = 0
                    self.median_cnv_list.clear()
                    for i in range(len(changepoint_list) - 1):
                        slice_start = (changepoint_list[i])  # offset
                        slice_end = (changepoint_list[i + 1])  # offset
                        slice_data = (combined_array[slice_start:slice_end])
                        if slice_data.size > 0:  # Proper way to check if numpy array is not empty
                            self.median_cnv_list.append(statistics.median(slice_data))
                        else:
                            # Handle the case for empty slice_data, e.g., append a default value or skip
                            self.median_cnv_list.append(0)  # Example: appending 0 or some other placeholder

                    # FOR CSV GENERATION
                    #
                    # for changepoint in changepoint_list:
                    #     self.changepoint_y.append(np_cnv[changepoint - 1 - cp_total])
                    # cp_total += len(np_cnv)
                    #
                    # for changepoint in changepoint_list:
                    #     self.changepoint_list_scaled.append(changepoint * self.cnv_dict["bin_width"])
                    # for element in changepoint_list:
                    #     self.contig_list.append(contig)
                    # self.changepoint_report.update({'chr': self.contig_list,'changepoint': self.changepoint_list_scaled,
                    #                                 'cnv': self.changepoint_y})
                    # print(len(self.contig_list))
                    # print(len(self.changepoint_y))
                    # print(len(self.changepoint_list_scaled))
                    # changepoint_report = pd.DataFrame(self.changepoint_report)
                    # changepoint_report.to_csv('/home/p2solo/cnsmeth-step_wise_startup/src/cnsmeth/subpages'
                    #                           '/changepoints.csv', index=False)
                    # self.download_cnv_data = 'changepoints.csv'

                    # Prepare data for plotting
                    data = list(
                        zip(
                            (np.arange(len(np_cnv)) + total) * self.cnv_dict["bin_width"],
                            cnv,
                        )
                    )
                    data = reduce_list(data)  # For faster processing?
                    total += len(cnv)  # Set start pos of ruptures coordinates to end of previous chromosome

                    print(f"running plots for {contig}")
                    # Update plot options, set min, max and series.
                    plot_to_update.options["xAxis"]["max"] = max
                    plot_to_update.options["xAxis"]["min"] = min
                    plot_to_update.options["series"].append(
                        {
                            "type": "scatter",
                            "name": contig,
                            "data": data,
                            "symbol": "circle" if self.scatter_toggle.value == "on" else "none",
                            "symbolSize": 5,
                            "markLine": {
                                "symbol": "none",
                                "data": [
                                    {
                                        "lineStyle": {"width": 2},
                                        "label": {"formatter": contig},
                                        "xAxis": (
                                                (total) * self.cnv_dict["bin_width"]

                                        ),
                                    },
                                ],
                            },
                        }
                    )

                    # Add new ruptures changepoint series

                    for element in changepoint_list:
                        if self.changepoint_toggle.value == "on":
                            plot_to_update.options["series"][0]["markLine"]["data"].append(
                                {
                                    "lineStyle": {"width": 1},
                                    "label": {"normal": {"show": False}},
                                    "name": contig,
                                    "xAxis": (
                                            (element) * self.cnv_dict["bin_width"]
                                    ),
                                }
                            )

                    if self.step_line_toggle.value == "on":
                        plot_to_update.options["series"].append(
                            {
                                "type": "line",
                                "name": f"Step Line {contig}",
                                "step": "end",
                                "data": list(
                                    zip([x * self.cnv_dict["bin_width"] for x in changepoint_list],
                                        self.median_cnv_list)),
                                "lineStyle": {"width": 2},

                            }
                        )
                    # print(list(zip([x * self.cnv_dict["bin_width"] for x in changepoint_list], self.median_cnv_list)))
                    # Add gene information to genevalueslist
                    for index, gene in self.gene_bed[
                        self.gene_bed["chrom"] == contig
                    ].iterrows():
                        genevalueslist[index] = f"{gene.chrom} - {gene.gene}"

            # If a specific chromosome is selected
            else:

                plot_to_update.options["series"] = []
                for counter, contig in enumerate(
                        natsort.natsorted(result.cnv), start=1
                ):
                    valueslist[counter] = contig
                # Retrieve CNV data from the chromosome selected based on the indexing of result.cnv
                contig, cnv = natsort.natsorted(result.cnv.items())[
                    int(self.chrom_filter) - 1
                    ]
                data = list(
                    zip((np.arange(len(cnv)) + total) * self.cnv_dict["bin_width"], cnv)
                )

                # Change point detection for the selected chromosome
                changepoint_list, changepoint_dict = self.detect_changepoints(contig, cnv, penalty, total,
                                                                              changepoint_dict)
                self.median_cnv_list.clear()
                for i in range(len(changepoint_list) - 1):
                    slice_start = changepoint_list[i]
                    slice_end = changepoint_list[i + 1]
                    slice_data = np.array(cnv[slice_start:slice_end])
                    if slice_data.size > 0:
                        self.median_cnv_list.append(statistics.median(slice_data))
                    else:
                        self.median_cnv_list.append(0)

                # Add final y-point to chrom, average of last 10 elements
                if np.array(cnv[-10:]).size > 0:
                    self.median_cnv_list.append(statistics.median(cnv[-10:]))
                else:
                    self.median_cnv_list.append(0)

                # Adjust min and max based on gene_target, no gene targets or gene_target for all chromosomes result in
                # Default min and max values (start to end of genome)
                if not gene_target:
                    min = min
                    max = "dataMax"

                else:
                    if gene_target == "All":
                        min = min
                        max = "dataMax"
                    else:
                        # Retrieve start and end positions of the gene
                        start_pos = self.gene_bed.iloc[int(gene_target)].start_pos
                        end_pos = self.gene_bed.iloc[int(gene_target)].end_pos
                        chrom = self.gene_bed.iloc[int(gene_target)].chrom
                        # counter is index
                        counter = 0
                        for counter, contig in enumerate(
                                natsort.natsorted(result.cnv), start=1
                        ):
                            # Checks if current chromosome matches selected chrom
                            if contig == chrom:
                                self.chrom_filter = counter
                                break  # Breaks out of loop once selected chromosome is found

                        min = start_pos - 10 * self.cnv_dict["bin_width"]
                        max = end_pos + 10 * self.cnv_dict["bin_width"]

                        # Ensure min and max are within certain restrictions
                        if start_pos - min > 2_000_000:
                            min = start_pos - 2_000_000
                        if max - end_pos > 2_000_000:
                            max = end_pos + 2_000_000
                        # Ensures min is not a negative value
                        if min < 0:
                            min = 0

                # Update plot
                plot_to_update.options["title"][
                    "text"
                ] = f"Copy Number Variation - {contig}"
                plot_to_update.options["xAxis"]["max"] = max
                plot_to_update.options["xAxis"]["min"] = min
                plot_to_update.options["series"].append(
                    {
                        "type": "scatter",
                        "name": contig,
                        "data": data,
                        "symbol": "circle" if self.scatter_toggle.value == "on" else "none",
                        "symbolSize": 10,
                        "itemStyle": {
                            "color": "rgb(136, 201, 253)"
                        },
                        "markLine": {
                            'symbol': ['none', 'none'],
                            "lineStyle": {
                                "color": "red",
                                "type": "dotted",

                            },
                            "data": [
                            ]
                        },
                        "markArea": {
                            "itemStyle": {"color": "rgba(255, 173, 177, 0.4)"},
                            "data": [],
                        },
                    }
                )

                # Add changepoints as markLine to the plot
                for changepoint in changepoint_dict.get(contig):
                    if self.changepoint_toggle.value == "on":
                        plot_to_update.options["series"][0]["markLine"]["data"].append(
                            {
                                "lineStyle": {"width": 3},
                                "label": {"normal": {"show": False}},
                                "name": contig,
                                "xAxis": (
                                        (changepoint) * self.cnv_dict["bin_width"]
                                ),
                            }
                        )

                # Add step line for the selected chromosome
                if self.step_line_toggle.value == "on":
                    plot_to_update.options["series"].append(
                        {
                            "type": "line",
                            "name": f"Step Line {contig}",
                            "step": "end",
                            "data": list(
                                zip([x * self.cnv_dict["bin_width"] for x in changepoint_list], self.median_cnv_list)),
                            "lineStyle": {
                                "width": 3,
                                "color": "black",
                            },

                        }
                    )

                for index, gene in self.gene_bed[
                    self.gene_bed["chrom"] == contig
                ].iterrows():
                    genevalueslist[index] = f"{gene.chrom} - {gene.gene}"

                for _, row in self.gene_bed[
                    self.gene_bed["chrom"] == contig
                ].iterrows():
                    plot_to_update.options["series"][0]["markArea"]["data"].append(
                        [
                            {
                                "name": row["gene"],
                                "xAxis": row["start_pos"],
                            },
                            {
                                "xAxis": row["end_pos"],
                            },
                        ]
                    )

            self.chrom_select.set_options(valueslist)
            self.gene_select.set_options(genevalueslist)
            plot_to_update.update()

    def show_previous_data(self, output):
        result = np.load(os.path.join(output, "CNV.npy"), allow_pickle="TRUE").item()
        self.result = Result(result)
        cnv_dict = np.load(
            os.path.join(output, "CNV_dict.npy"), allow_pickle=True
        ).item()
        self.cnv_dict["bin_width"] = cnv_dict["bin_width"]
        self.cnv_dict["variance"] = cnv_dict["variance"]
        self._update_cnv_plot()


def test_me(
        port: int,
        threads: int,
        watchfolder: str,
        output: str,
        reload: bool = False,
        browse: bool = False,
        target_panel: str = "rCNS2",
):
    my_connection = None
    app.add_static_files("/fonts", str(Path(__file__).parent / "../fonts"))
    with theme.frame("Copy Number Variation Testing."):
        TestObject = CNVAnalysis(
            threads,
            output,
            progress=True,
            # bamqueue=self.bamforcnv,
            # summary=cnvsummary,
            target_panel=target_panel,
        )
        # TestObject = CNVAnalysis(threads, output, progress=True)
    if not browse:
        path = watchfolder
        searchdirectory = os.fsencode(path)
        for root, d_names, f_names in os.walk(searchdirectory):
            directory = os.fsdecode(root)
            for f in f_names:
                filename = os.fsdecode(f)
                if filename.endswith(".bam"):
                    TestObject.add_bam(os.path.join(directory, filename))
    else:
        # print("Browse mode not implemented.")
        TestObject.progress_trackers.visible = False
        TestObject.show_previous_data(output)
    ui.run(port=port, reload=reload)


@click.command()
@click.option(
    "--port",
    default=22345,
    help="Port for GUI",
)
@click.option("--threads", default=4, help="Number of threads available.")
@click.argument(
    "watchfolder",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.argument(
    "output",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.option(
    "--browse",
    is_flag=True,
    show_default=True,
    default=False,
    help="Browse Historic Data.",
)
@click.option(
    "--target_panel",
    "-t",
    default="rCNS2",
    help="Select analysis gene panel from one of these options. Default is rCNS2",
    type=click.Choice(
        ["rCNS2", "AML"],
        case_sensitive=True,
    ),
)
def main(port, threads, watchfolder, output, browse, target_panel):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    if browse:
        # Handle the case when --browse is set
        click.echo("Browse mode is enabled. Only the output folder is required.")
        test_me(
            port=port,
            reload=False,
            threads=threads,
            # simtime=simtime,
            watchfolder=None,
            output=watchfolder,
            # sequencing_summary=sequencing_summary,
            # showerrors=showerrors,
            browse=browse,
            # exclude=exclude,
        )
        # Your logic for browse mode
    else:
        # Handle the case when --browse is not set
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if watchfolder is None or output is None:
            click.echo("Watchfolder and output are required when --browse is not set.")
            sys.exit(1)
        test_me(
            port=port,
            reload=False,
            threads=threads,
            # simtime=simtime,
            watchfolder=watchfolder,
            output=output,
            # sequencing_summary=sequencing_summary,
            # showerrors=showerrors,
            browse=browse,
            target_panel=target_panel,
            # exclude=exclude,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    main()

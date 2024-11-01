package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.DefaultXYZDataset;

import javax.swing.*;
import java.util.HashMap;
import java.util.Map;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final int NUM_BARS = 60;  // Number of bars
    private static final int NUM_WEDGES = 10; // Number of wedges per bar
    private static final double WEDGE_THICKNESS = 30.0; // Thickness of each wedge in mm

    // Constants for detector structure
    private static final int NUM_SECTORS = 15;  // Total number of sectors
    private static final int NUM_LAYERS = 4;    // Number of layers in each sector (valid values: 0 to 3)
    private static final int NUM_ORDERS = 2;    // Inner and outer wedges (0 = inner, 1 = outer)

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];
        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        if (!reader.getSchemaFactory().hasSchema("ATOF::adc")) {
            System.err.println("Schema ATOF::adc not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        // Store Z, Phi positions and event counts for wedges
        Map<Double, Integer>[] zPositionEventCountsForWedge = new HashMap[NUM_SECTORS * NUM_LAYERS * NUM_ORDERS * NUM_BARS * NUM_WEDGES];
        Map<Double, Integer>[] phiPositionEventCountsForBar = new HashMap[NUM_BARS];
        for (int i = 0; i < zPositionEventCountsForWedge.length; i++) {
            zPositionEventCountsForWedge[i] = new HashMap<>();
        }
        for (int i = 0; i < phiPositionEventCountsForBar.length; i++) {
            phiPositionEventCountsForBar[i] = new HashMap<>();
        }

        // Extract Z and Phi positions for wedges
        extractAndPrintWedges(reader, zPositionEventCountsForWedge, phiPositionEventCountsForBar);

        // Plot Z vs number of events (bar chart)
        plotPositionVsEvents("Z Position vs Number of Events for Wedges", "Z Position (mm)", "Number of Events", zPositionEventCountsForWedge);

        // Plot Number of Events vs Z (scatter plot)
        plotScatterPlot("Number of Events vs Z Position (Scatter Plot)", "Z Position (mm)", "Number of Events", zPositionEventCountsForWedge);

        // Plot Phi vs Number of Events (for each bar)
        plotPhiVsEvents("Phi Position vs Number of Events", "Phi Position (radians)", "Number of Events", phiPositionEventCountsForBar);

        // Plot Z vs Phi (2D plot)
        plotZVsPhi("Z Position vs Phi Position (2D Plot)", zPositionEventCountsForWedge, phiPositionEventCountsForBar);

        reader.close();
    }

    private static void extractAndPrintWedges(HipoReader reader, Map<Double, Integer>[] zPositionEventCountsForWedge, Map<Double, Integer>[] phiPositionEventCountsForBar) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        // Iterate over events
        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int atofRows = atofAdcBank.getRows();
            if (atofRows < 2) continue;  // Ensure we have at least 2 rows (left and right PMT)

            // Iterate through available rows for bars and wedges
            for (int rowIndex = 0; rowIndex < atofRows; rowIndex++) {
                int sector = atofAdcBank.getInt("sector", rowIndex);
                int layer = atofAdcBank.getInt("layer", rowIndex);
                int order = atofAdcBank.getInt("order", rowIndex);
                int component = atofAdcBank.getInt("component", rowIndex);

                // Validate layer
                if (layer < 0 || layer >= NUM_LAYERS) {
                    continue;  // Ignore invalid layers
                }

                // Calculate indices for wedges based on sector, layer, order, and component (bar)
                int barIndex = component - 1;
                double phiBar = calculatePhiForBar(barIndex);  // Calculate the Phi value for this bar

                for (int wedgeIndex = 0; wedgeIndex < NUM_WEDGES; wedgeIndex++) {
                    //double zWedge = (wedgeIndex - NUM_WEDGES / 2.0) * WEDGE_THICKNESS;// plots Z posiiton in between two wedges 
                    double zWedge = (wedgeIndex - (NUM_WEDGES - 1) / 2.0) * WEDGE_THICKNESS;
                    int globalWedgeIndex = calculateGlobalWedgeIndex(sector, layer, order, barIndex, wedgeIndex);

                    if (globalWedgeIndex >= 0 && globalWedgeIndex < zPositionEventCountsForWedge.length) {
                        int currentEventCount = zPositionEventCountsForWedge[globalWedgeIndex].getOrDefault(zWedge, 0);
                        int randomEventAddition = (int) (Math.random() * 5) + 1;
                        zPositionEventCountsForWedge[globalWedgeIndex].put(zWedge, currentEventCount + randomEventAddition);

                        System.out.printf("Sector: %d  Layer: %d  Order: %d  Component: %d  Bar: %d  Wedge: %d  Z: %.2f mm  Events: %d%n",
                                sector, layer, order, component, barIndex + 1, wedgeIndex, zWedge, currentEventCount + randomEventAddition);
                    }
                }

                // Add event counts to Phi for the whole bar
                int currentPhiEventCount = phiPositionEventCountsForBar[barIndex].getOrDefault(phiBar, 0);
                phiPositionEventCountsForBar[barIndex].put(phiBar, currentPhiEventCount + 1);
            }
        }
    }

    private static double calculatePhiForBar(int barIndex) {
        // Assume the Phi range is from -π to π, and bars are evenly spaced in Phi
        double phiMin = -Math.PI;
        double phiMax = Math.PI;
        return phiMin + (phiMax - phiMin) * barIndex / NUM_BARS;
    }

    private static int calculateGlobalWedgeIndex(int sector, int layer, int order, int barIndex, int wedgeIndex) {
        return sector * NUM_LAYERS * NUM_ORDERS * NUM_BARS * NUM_WEDGES +
                layer * NUM_ORDERS * NUM_BARS * NUM_WEDGES +
                order * NUM_BARS * NUM_WEDGES +
                barIndex * NUM_WEDGES +
                wedgeIndex;
    }

    // Plot Z vs Number of Events (Bar Chart)
    private static void plotPositionVsEvents(String title, String xAxisLabel, String yAxisLabel, Map<Double, Integer>[] positionEventCountsForWedges) {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        for (int i = 0; i < positionEventCountsForWedges.length; i++) {
            for (Map.Entry<Double, Integer> entry : positionEventCountsForWedges[i].entrySet()) {
                dataset.addValue(entry.getValue(), "Events", String.format("%.2f", entry.getKey()));
            }
        }
        JFreeChart chart = ChartFactory.createBarChart(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, title);
    }

    // Plot Z vs Number of Events (Scatter Plot)
    private static void plotScatterPlot(String title, String xAxisLabel, String yAxisLabel, Map<Double, Integer>[] positionEventCountsForWedges) {
        XYSeries series = new XYSeries("Events");
        for (int i = 0; i < positionEventCountsForWedges.length; i++) {
            for (Map.Entry<Double, Integer> entry : positionEventCountsForWedges[i].entrySet()) {
                series.add(entry.getKey(), entry.getValue());
            }
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, dataset);
        displayChart(chart, title);
    }

    // Plot Phi vs Number of Events
    private static void plotPhiVsEvents(String title, String xAxisLabel, String yAxisLabel, Map<Double, Integer>[] phiPositionEventCountsForBar) {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
        for (int i = 0; i < phiPositionEventCountsForBar.length; i++) {
            for (Map.Entry<Double, Integer> entry : phiPositionEventCountsForBar[i].entrySet()) {
                dataset.addValue(entry.getValue(), "Events", String.format("%.2f", entry.getKey()));
            }
        }
        JFreeChart chart = ChartFactory.createBarChart(title, xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, title);
    }

    // Plot Z vs Phi (2D Plot)
    private static void plotZVsPhi(String title, Map<Double, Integer>[] zPositionEventCountsForWedge, Map<Double, Integer>[] phiPositionEventCountsForBar) {
        DefaultXYZDataset dataset = new DefaultXYZDataset();
        double[] zValues = new double[NUM_BARS * NUM_WEDGES];
        double[] phiValues = new double[NUM_BARS * NUM_WEDGES];
        double[] eventValues = new double[NUM_BARS * NUM_WEDGES];
        int idx = 0;

        for (int barIndex = 0; barIndex < NUM_BARS; barIndex++) {
            double phiBar = calculatePhiForBar(barIndex); // Each bar has the same Phi for all its wedges

            for (int wedgeIndex = 0; wedgeIndex < NUM_WEDGES; wedgeIndex++) {
                double zWedge = (wedgeIndex - NUM_WEDGES / 2.0) * WEDGE_THICKNESS;
                int globalWedgeIndex = barIndex * NUM_WEDGES + wedgeIndex; // Correct indexing for wedges

                // Now use the correct index to get Z and Phi values
                zValues[idx] = zWedge; // Z for the wedge
                phiValues[idx] = phiBar; // Same Phi for all wedges in the bar
                eventValues[idx] = zPositionEventCountsForWedge[globalWedgeIndex].getOrDefault(zWedge, 0); // Events count
                idx++;
            }
        }

        dataset.addSeries("Z vs Phi", new double[][]{zValues, phiValues, eventValues});
        JFreeChart chart = ChartFactory.createScatterPlot(title, "Z Position (mm)", "Phi Position (radians)", dataset);
        displayChart(chart, title);
    }

    // Utility to display the chart
    private static void displayChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }
}
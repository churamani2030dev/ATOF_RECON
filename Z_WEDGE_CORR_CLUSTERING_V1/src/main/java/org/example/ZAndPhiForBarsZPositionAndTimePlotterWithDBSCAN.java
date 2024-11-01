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

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final int NUM_BARS = 60;  // Number of bars
    private static final int NUM_WEDGES = 10; // Number of wedges per bar
    private static final double WEDGE_SPACING = 30.0; // Wedge spacing (mm)
    private static final double VELOCITY_EFF = 200.0;
    private static final double Z_THRESHOLD = 5.0;
    private static final double PHI_THRESHOLD = 0.01;
    private static final double TIME_THRESHOLD = 1.0;

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

        processEvents(reader);
        reader.close();
    }

    private static void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();

        List<Integer> clusterSizes = new ArrayList<>();
        List<Integer> clusterIndices = new ArrayList<>();

        List<Double> deltaZList = new ArrayList<>();
        List<Double> deltaPhiList = new ArrayList<>();
        List<Double> deltaTimeList = new ArrayList<>();
        List<Integer> eventIndicesForDeltas = new ArrayList<>();

        List<Double> zBarValues = new ArrayList<>();
        List<Double> zWedgeValues = new ArrayList<>();
        List<Double> phiBarValues = new ArrayList<>();
        List<Double> timeBarValues = new ArrayList<>();
        List<Double> timeWedgeValues = new ArrayList<>();
        List<Integer> eventIndices = new ArrayList<>();

        int eventIndex = 0;
        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            int numRows = atofAdcBank.getRows();
            System.out.println("\nProcessing a new event with " + numRows + " hits...");

            List<EventData> eventsData = new ArrayList<>();
            List<Cluster> clusters = new ArrayList<>();

            for (int hitIndex = 0; hitIndex < numRows; hitIndex++) {
                int sector = atofAdcBank.getInt("sector", hitIndex);
                int layer = atofAdcBank.getInt("layer", hitIndex);
                int component = atofAdcBank.getInt("component", hitIndex);
                int order = atofAdcBank.getInt("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                short pedestal = atofAdcBank.getShort("ped", hitIndex);
                double time = atofAdcBank.getFloat("time", hitIndex);

                double wedgeZ = calculateZForWedge(hitIndex);
                double barZ = calculateZForBar(atofAdcBank, hitIndex);
                double phi = calculatePhiForBar(hitIndex);
                double barTime = calculateBarTime(atofAdcBank, hitIndex);

                eventsData.add(new EventData(wedgeZ, barZ, phi, phi, time, barTime, sector, layer, component, order, adc, pedestal));

                // Print detailed hit information
                System.out.printf("Hit %d -> Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Pedestal: %d, ZWedge: %.2f, ZBar: %.2f, Phi: %.2f, TimeW: %.2f, TimeB: %.2f\n",
                        hitIndex, sector, layer, component, order, adc, pedestal, wedgeZ, barZ, phi, time, barTime);

                // Add deltas for plotting
                deltaZList.add(Math.abs(barZ - wedgeZ));
                deltaPhiList.add(Math.abs(phi - phi)); // Phi difference (Bar - Wedge)
                deltaTimeList.add(Math.abs(barTime - time)); // Time difference (Bar - Wedge)

                zBarValues.add(barZ);
                zWedgeValues.add(wedgeZ);
                phiBarValues.add(phi);
                timeBarValues.add(barTime);
                timeWedgeValues.add(time);
                eventIndices.add(eventIndex);
            }

            formClusters(eventsData, clusters);

            if (!clusters.isEmpty()) {
                System.out.println("Clusters formed in this event:");
                printClusters(clusters);

                clusterSizes.add(clusters.size());
                clusterIndices.add(eventIndex);
            }

            eventIndicesForDeltas.add(eventIndex);
            eventIndex++;
        }

        // Plot clusters, deltas, and event-based data
        plotClusterSizeVsEventIndex(clusterSizes, clusterIndices);
        plotDeltaVsEventIndex(deltaZList, eventIndicesForDeltas, "Delta Z vs Event Index", "Delta Z");
        plotDeltaVsEventIndex(deltaPhiList, eventIndicesForDeltas, "Delta Phi vs Event Index", "Delta Phi");
        plotDeltaVsEventIndex(deltaTimeList, eventIndicesForDeltas, "Delta Time vs Event Index", "Delta Time");

        plotZBarVsEventIndex(zBarValues, eventIndices);
        plotZWedgeVsEventIndex(zWedgeValues, eventIndices);
        plotPhiBarVsEventIndex(phiBarValues, eventIndices);
        plotTimeVsEventIndex(timeBarValues, eventIndices, "TimeBar");
        plotTimeVsEventIndex(timeWedgeValues, eventIndices, "TimeWedge");
    }

    // Method to plot Cluster Size vs Event Index
    private static void plotClusterSizeVsEventIndex(List<Integer> clusterSizes, List<Integer> clusterIndices) {
        XYSeries series = new XYSeries("Cluster Size vs Event Index");
        for (int i = 0; i < clusterSizes.size(); i++) {
            series.add(clusterIndices.get(i), clusterSizes.get(i));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                "Cluster Size vs Event Index", "Event Index", "Cluster Size",
                dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, "Cluster Size vs Event Index Plot");
    }

    // Method to plot Delta values vs Event Index
    private static void plotDeltaVsEventIndex(List<Double> deltaValues, List<Integer> eventIndices, String plotTitle, String yAxisLabel) {
        XYSeries series = new XYSeries(plotTitle);
        for (int i = 0; i < deltaValues.size(); i++) {
            series.add(eventIndices.get(i), deltaValues.get(i));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                plotTitle, "Event Index", yAxisLabel,
                dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, plotTitle);
    }

    // Method to plot ZBar vs Event Index
    private static void plotZBarVsEventIndex(List<Double> ZBarValues, List<Integer> eventIndices) {
        XYSeries series = new XYSeries("ZBar vs Event Index");
        for (int i = 0; i < ZBarValues.size(); i++) {
            series.add(eventIndices.get(i), ZBarValues.get(i));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                "ZBar vs Event Index", "Event Index", "ZBar (mm)",
                dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, "ZBar vs Event Index Plot");
    }

    // Method to plot ZWedge vs Event Index
    private static void plotZWedgeVsEventIndex(List<Double> ZWedgeValues, List<Integer> eventIndices) {
        XYSeries series = new XYSeries("ZWedge vs Event Index");
        for (int i = 0; i < ZWedgeValues.size(); i++) {
            series.add(eventIndices.get(i), ZWedgeValues.get(i));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                "ZWedge vs Event Index", "Event Index", "ZWedge (mm)",
                dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, "ZWedge vs Event Index Plot");
    }


		         // Method to plot PhiBar vs Event Index
    private static void plotPhiBarVsEventIndex(List<Double> PhiBarValues, List<Integer> eventIndices) {
        XYSeries series = new XYSeries("PhiBar vs Event Index");
        for (int i = 0; i < PhiBarValues.size(); i++) {
            series.add(eventIndices.get(i), PhiBarValues.get(i));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                "PhiBar vs Event Index", "Event Index", "PhiBar (radians)",
                dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, "PhiBar vs Event Index Plot");
    }

    // Method to plot Time values (Bar or Wedge) vs Event Index
    private static void plotTimeVsEventIndex(List<Double> TimeValues, List<Integer> eventIndices, String timeType) {
        XYSeries series = new XYSeries(timeType + " vs Event Index");
        for (int i = 0; i < TimeValues.size(); i++) {
            series.add(eventIndices.get(i), TimeValues.get(i));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                timeType + " vs Event Index", "Event Index", timeType + " (ns)",
                dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, timeType + " vs Event Index Plot");
    }

    // Utility to display the chart
    private static void displayChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }

    // EventData class to store information about hits
    private static class EventData {
        double zWedge, zBar, phiWedge, phiBar, timeWedge, timeBar;
        int sector, layer, component, order, adc;
        short pedestal;

        EventData(double zWedge, double zBar, double phiWedge, double phiBar,
                  double timeWedge, double timeBar, int sector, int layer,
                  int component, int order, int adc, short pedestal) {
            this.zWedge = zWedge;
            this.zBar = zBar;
            this.phiWedge = wrapPhi(phiWedge);
            this.phiBar = wrapPhi(phiBar);
            this.timeWedge = timeWedge;
            this.timeBar = timeBar;
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.pedestal = pedestal;
        }
    }

    // Cluster class to handle event clustering
    private static class Cluster {
        List<EventData> events = new ArrayList<>();

        public void addEvent(EventData event) {
            events.add(event);
        }

        public int getClusterSize() {
            return events.size();
        }

        public boolean isValidCluster() {
            return events.size() >= 2; // At least 2 hits for a valid cluster
        }

        public void printCluster(int clusterIndex) {
            System.out.printf("Cluster %d - Size: %d\n", clusterIndex, getClusterSize());
            for (EventData event : events) {
                System.out.printf("  ZW: %.2f, ZB: %.2f, PhiW: %.2f, PhiB: %.2f, TW: %.2f, TB: %.2f, Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Pedestal: %d\n",
                        event.zWedge, event.zBar, event.phiWedge, event.phiBar, event.timeWedge, event.timeBar,
                        event.sector, event.layer, event.component, event.order, event.adc, event.pedestal);
            }
        }
    }

    // Method to form clusters based on proximity criteria
    private static void formClusters(List<EventData> eventsData, List<Cluster> clusters) {
        for (EventData event : eventsData) {
            boolean addedToCluster = false;

            for (Cluster cluster : clusters) {
                if (isWithinProximity(cluster, event)) {
                    cluster.addEvent(event);
                    addedToCluster = true;
                    break;
                }
            }

            if (!addedToCluster) {
                Cluster newCluster = new Cluster();
                newCluster.addEvent(event);
                clusters.add(newCluster);
            }
        }
    }

    // Method to determine if an event belongs to an existing cluster
    private static boolean isWithinProximity(Cluster cluster, EventData event) {
        for (EventData clusterEvent : cluster.events) {
            if (Math.abs(clusterEvent.zBar - event.zWedge) < Z_THRESHOLD &&
                Math.abs(clusterEvent.phiBar - event.phiWedge) < PHI_THRESHOLD &&
                Math.abs(clusterEvent.timeBar - event.timeWedge) < TIME_THRESHOLD) {
                return true;
            }
        }
        return false;
    }

    // Method to print clusters
    private static void printClusters(List<Cluster> clusters) {
        int clusterIndex = 1;
        for (Cluster cluster : clusters) {
            if (cluster.isValidCluster()) {
                cluster.printCluster(clusterIndex++);
            }
        }
    }

    // Utility methods for calculating Z and Phi values
    private static double calculateZForBar(Bank bank, int rowIndex) {
        double timeLeftPMT = bank.getFloat("time", rowIndex);
        double timeRightPMT = bank.getFloat("time", rowIndex + 1);
        return VELOCITY_EFF * (timeRightPMT - timeLeftPMT) / 2.0;
    }

    private static double calculateBarTime(Bank bank, int rowIndex) {
        double timeLeftPMT = bank.getFloat("time", rowIndex);
        double timeRightPMT = bank.getFloat("time", rowIndex + 1);
        return Math.min(timeLeftPMT, timeRightPMT);
    }

    private static double calculateZForWedge(int wedgeIndex) {
        return (wedgeIndex - (NUM_WEDGES - 1) / 2.0) * WEDGE_SPACING;
    }

    private static double calculatePhiForBar(int barIndex) {
        double phi = -Math.PI + (2 * Math.PI) * barIndex / NUM_BARS;
        return wrapPhi(phi);
    }

    // Method to wrap phi values between -PI and +PI
    private static double wrapPhi(double phi) {
        while (phi > Math.PI) phi -= 2 * Math.PI;
        while (phi < -Math.PI) phi += 2 * Math.PI;
        return phi;
    }
}

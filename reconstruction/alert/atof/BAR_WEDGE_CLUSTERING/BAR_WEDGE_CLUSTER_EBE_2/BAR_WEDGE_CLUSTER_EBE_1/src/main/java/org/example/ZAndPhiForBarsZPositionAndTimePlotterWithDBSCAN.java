package org.example;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final int NUM_WEDGES = 10;
    private static final int NUM_BARS = 60;
    private static final int MAX_LAYER = 3;
    private static final double WEDGE_SPACING = 30.0;
    private static final double VELOCITY_EFF = 200.0;
    private static final double Z_THRESHOLD = 5.0;
    private static final double PHI_THRESHOLD = 0.01;
    private static final double TIME_THRESHOLD = 1.0;

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

    private static class Cluster {
        List<EventData> events = new ArrayList<>();

        public void addEvent(EventData event) {
            events.add(event);
        }

        public int getClusterSize() {
            return events.size();
        }

        public boolean isValidCluster() {
            return events.size() >= 2;
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

        // Loop over events and perform clustering per event
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

                if (layer < 0 || layer > MAX_LAYER) {
                    System.out.println("Invalid layer: " + layer + ". Skipping hit.");
                    continue;
                }

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

                System.out.printf("Hit %d -> ZW: %.2f, ZB: %.2f, Phi: %.2f, TimeW: %.2f, TimeB: %.2f, Sector: %d, Layer: %d, Component: %d\n",
                        hitIndex, wedgeZ, barZ, phi, time, barTime, sector, layer, component);

                // Add deltas for plotting
                deltaZList.add(Math.abs(barZ - wedgeZ));
                deltaPhiList.add(Math.abs(phi - phi));  // No delta between phiBar and phiWedge since they're the same
                deltaTimeList.add(Math.abs(barTime - time));
                eventIndicesForDeltas.add(eventIndex);
            }

            formClusters(eventsData, clusters);

            if (!clusters.isEmpty()) {
                System.out.println("Clusters formed in this event:");
                printClusters(clusters);

                clusterSizes.add(clusters.size());
                clusterIndices.add(eventIndex);
            }

            eventIndex++;
        }

        // Plotting the results for cluster size and deltas
        plotClusterSizeVsEventIndex(clusterSizes, clusterIndices);
        plotDeltaVsEventIndex(deltaZList, eventIndicesForDeltas, "Delta Z vs Event Index", "Delta Z");
        plotDeltaVsEventIndex(deltaPhiList, eventIndicesForDeltas, "Delta Phi vs Event Index", "Delta Phi");
        plotDeltaVsEventIndex(deltaTimeList, eventIndicesForDeltas, "Delta Time vs Event Index", "Delta Time");
    }

    private static void plotClusterSizeVsEventIndex(List<Integer> clusterSizes, List<Integer> clusterIndices) {
        XYSeries series = new XYSeries("Cluster Size vs Event Index");
        for (int i = 0; i < clusterSizes.size(); i++) {
            series.add(clusterIndices.get(i), clusterSizes.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                "Cluster Size vs Event Index",
                "Event Index", "Cluster Size",
                dataset, PlotOrientation.VERTICAL, true, true, false);

        displayChart(chart, "Cluster Size vs Event Index Plot");
    }

    private static void plotDeltaVsEventIndex(List<Double> deltaValues, List<Integer> eventIndices, String plotTitle, String yAxisLabel) {
        XYSeries series = new XYSeries(plotTitle);
        for (int i = 0; i < deltaValues.size(); i++) {
            series.add(eventIndices.get(i), deltaValues.get(i));
        }

        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(
                plotTitle,
                "Event Index", yAxisLabel,
                dataset, PlotOrientation.VERTICAL, true, true, false);

        displayChart(chart, plotTitle);
    }

    private static void displayChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }

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

    private static double wrapPhi(double phi) {
        while (phi > Math.PI) phi -= 2 * Math.PI;
        while (phi < -Math.PI) phi += 2 * Math.PI;
        return phi;
    }

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

    private static void printClusters(List<Cluster> clusters) {
        int clusterIndex = 1;
        for (Cluster cluster : clusters) {
            if (cluster.isValidCluster()) {
                cluster.printCluster(clusterIndex++);
            }
        }
    }
}

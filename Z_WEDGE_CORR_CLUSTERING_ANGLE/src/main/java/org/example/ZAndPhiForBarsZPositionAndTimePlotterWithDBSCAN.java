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

    private static final int NUM_WEDGES = 10;  // 10 wedges per bar
    private static final int NUM_BARS = 60;    // 60 bars total
    private static final double WEDGE_SPACING = 30.0;  // Wedge spacing in mm
    private static final double VELOCITY_EFF = 200.0;  // Effective velocity for time calculations

    private static class EventData {
        Double zWedge = null, zBar = null, phiWedge = null, phiBar = null, timeWedge = null, timeBar = null;
        int sector, layer, component, order, adc;
        short pedestal;

        EventData(Double zWedge, Double zBar, Double phiWedge, Double phiBar,
                  Double timeWedge, Double timeBar, int sector, int layer,
                  int component, int order, int adc, short pedestal) {
            this.zWedge = zWedge;
            this.zBar = zBar;
            this.phiWedge = phiWedge;
            this.phiBar = phiBar;
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
        EventData wedgeHit = null;
        EventData barHitLeft = null;
        EventData barHitRight = null;

        public void addEvent(EventData event) {
            if (event.zBar != null && event.order == 0) {
                barHitLeft = event;
            } else if (event.zBar != null && event.order == 1) {
                barHitRight = event;
            } else if (event.zWedge != null) {
                wedgeHit = event;
            }
            events.add(event);
        }

        public boolean isValidCluster() {
            return (barHitLeft != null && barHitRight != null && wedgeHit != null);
        }

        public void printCluster(int clusterIndex) {
            System.out.printf("Cluster %d - Size: %d\n", clusterIndex, events.size());
            for (EventData event : events) {
                System.out.printf("  %s %s %s %s TW: %.2f, TB: %.2f, Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Pedestal: %d\n",
                        (event.zWedge != null ? String.format("ZW: %.2f", event.zWedge) : "ZW: N/A"),
                        (event.zBar != null ? String.format("ZB: %.2f", event.zBar) : "ZB: N/A"),
                        (event.phiWedge != null ? String.format("PhiW: %.2f", event.phiWedge) : "PhiW: N/A"),
                        (event.phiBar != null ? String.format("PhiB: %.2f", event.phiBar) : "PhiB: N/A"),
                        (event.timeWedge != null ? event.timeWedge : 0.0),
                        (event.timeBar != null ? event.timeBar : 0.0),
                        event.sector, event.layer, event.component, event.order, event.adc, event.pedestal);
            }

            // Calculate delta values
            if (isValidCluster()) {
                double deltaZ = barHitLeft.zBar - wedgeHit.zWedge;  // Z of the bar vs Z of the wedge
                double deltaPhi = barHitLeft.phiBar - wedgeHit.phiWedge;  // Phi of the bar vs Phi of the wedge
                double deltaTime = barHitLeft.timeBar - wedgeHit.timeWedge;  // Time of the bar vs Time of the wedge
                System.out.printf("Delta Z: %.2f, Delta Phi: %.2f, Delta Time: %.2f\n", deltaZ, deltaPhi, deltaTime);
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
            Cluster currentCluster = new Cluster();

            for (int hitIndex = 0; hitIndex < numRows; hitIndex++) {
                int sector = atofAdcBank.getInt("sector", hitIndex);
                int layer = atofAdcBank.getInt("layer", hitIndex);
                int component = atofAdcBank.getInt("component", hitIndex);
                int order = atofAdcBank.getInt("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                short pedestal = atofAdcBank.getShort("ped", hitIndex);
                double time = atofAdcBank.getFloat("time", hitIndex);

                // Handle null layer identification and default to layer 0 for invalid values
                if (layer >= 16) {
                    layer = layer == 16 ? 2 : 3;
                }

                Double zWedge = null, zBar = null, phi = null, barTime = null, wedgeTime = null;

                if (layer == 0) {
                    zBar = Double.valueOf(calculateZForBar(atofAdcBank, hitIndex));
                    phi = calculatePhiForBar(component - 1);
                    barTime = Double.valueOf(atofAdcBank.getFloat("time", hitIndex));
                } else if (layer == 2 || layer == 3) {
                    zWedge = Double.valueOf(calculateZForWedge(component));
                    phi = calculatePhiForBar(component - 1);
                    wedgeTime = Double.valueOf(atofAdcBank.getFloat("time", hitIndex));
                }

                currentCluster.addEvent(new EventData(zWedge, zBar, phi, phi, wedgeTime, barTime, sector, layer, component, order, adc, pedestal));
            }

            if (currentCluster.isValidCluster()) {
                clusters.add(currentCluster);
                System.out.println("Clusters formed in this event:");
                currentCluster.printCluster(clusters.size());

                deltaZList.add(currentCluster.barHitLeft.zBar - currentCluster.wedgeHit.zWedge);
                deltaPhiList.add(currentCluster.barHitLeft.phiBar - currentCluster.wedgeHit.phiWedge);
                deltaTimeList.add(currentCluster.barHitLeft.timeBar - currentCluster.wedgeHit.timeWedge);
                eventIndicesForDeltas.add(eventIndex);
            }

            if (!clusters.isEmpty()) {
                clusterSizes.add(clusters.size());
                clusterIndices.add(eventIndex);
            }

            eventIndex++;
        }

        // Plot cluster sizes
        plotClusterSizeVsEventIndex(clusterSizes, clusterIndices);

        // Plot delta Z, delta Phi, delta Time
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

    // Calculation methods for Z and Phi values
    private static double calculateZForBar(Bank bank, int rowIndex) {
        double timeLeftPMT = bank.getFloat("time", rowIndex);
        double timeRightPMT = bank.getFloat("time", rowIndex + 1);
        return VELOCITY_EFF * (timeRightPMT - timeLeftPMT) / 2.0;
    }

    private static double calculateZForWedge(int wedgeIndex) {
        return (wedgeIndex - (NUM_WEDGES - 1) / 2.0) * WEDGE_SPACING;
    }

    private static double calculatePhiForBar(int barIndex) {
        double phiMin = -Math.PI;
        double phiMax = Math.PI;
        return phiMin + (phiMax - phiMin) * barIndex / NUM_BARS;
    }
}

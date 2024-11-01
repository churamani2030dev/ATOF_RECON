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

public class ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN{

    private static final double WEDGE_SPACING = 30.0;
    private static final double VELOCITY_EFF = 200.0; // Effective velocity in cm/ns
    private static final int NUM_BARS = 60;
    private static final int NUM_WEDGES = 600;

    private static final double DELTA_Z_THRESHOLD = 5.0;
    private static final double DELTA_TIME_THRESHOLD = 1.0;
    private static final double DELTA_PHI_THRESHOLD = 0.1;

    private List<Double> deltaZList = new ArrayList<>();
    private List<Double> deltaTimeList = new ArrayList<>();
    private List<Double> deltaPhiList = new ArrayList<>();
    private List<Integer> clusterSizes = new ArrayList<>();
    private List<Integer> eventIndices = new ArrayList<>();

    private static class EventData {
        Double z, phi;
        double time;
        int sector, layer, component, order, adc;
        short pedestal;

        EventData(Double z, Double phi, double time, int sector, int layer, int component, int order, int adc, short pedestal) {
            this.z = z;
            this.phi = phi;
            this.time = time;
            this.sector = sector;
            this.layer = layer;
            this.component = component;
            this.order = order;
            this.adc = adc;
            this.pedestal = pedestal;
        }

        public void printDetails(String type) {
            System.out.printf("%s Hit -> Z: %.2f, Phi: %.2f, Time: %.2f, Sector: %d, Layer: %d, Component: %d, Order: %d, ADC: %d, Pedestal: %d\n",
                    type, z != null ? z : 0.0, phi != null ? phi : 0.0, time, sector, layer, component, order, adc, pedestal);
        }
    }

    private class Cluster {
        EventData wedgeHit = null;
        EventData barHitLeft = null;
        EventData barHitRight = null;

        public void addHit(EventData hit) {
            if (hit.layer == 2 || hit.layer == 3) {
                wedgeHit = hit;
            } else if (hit.layer == 0 && hit.order == 0) {
                barHitLeft = hit;
            } else if (hit.layer == 0 && hit.order == 1) {
                barHitRight = hit;
            }
        }

        public boolean isQualifiedCluster() {
            return wedgeHit != null && barHitLeft != null && barHitRight != null;
        }

        public void printAndRecordCluster(int clusterIndex) {
            if (!isQualifiedCluster()) {
                System.out.printf("Cluster %d is incomplete and cannot be processed (missing bar or wedge hits).\n", clusterIndex);
                return;
            }

            if (barHitLeft.z == null || barHitRight.z == null || wedgeHit.z == null) {
                System.out.printf("Cluster %d has null values in Z positions and cannot be processed.\n", clusterIndex);
                return;
            }
            if (barHitLeft.phi == null || wedgeHit.phi == null) {
                System.out.printf("Cluster %d has null values in Phi positions and cannot be processed.\n", clusterIndex);
                return;
            }

            double zBar = calculateZForBar(barHitLeft.time, barHitRight.time);
            double phiBar = barHitLeft.phi;
            double timeBar = Math.min(barHitLeft.time, barHitRight.time);

            double deltaZ = Math.abs(zBar - wedgeHit.z);
            double deltaPhi = Math.abs(phiBar - wedgeHit.phi);
            double deltaTime = Math.abs(timeBar - wedgeHit.time);

            System.out.printf("Cluster %d:\n", clusterIndex);
            System.out.printf("  Bar Avg Z: %.2f, Wedge Z: %.2f, Delta Z: %.2f\n", zBar, wedgeHit.z, deltaZ);
            System.out.printf("  Bar Phi: %.2f, Wedge Phi: %.2f, Delta Phi: %.2f\n", phiBar, wedgeHit.phi, deltaPhi);
            System.out.printf("  Bar Time: %.2f, Wedge Time: %.2f, Delta Time: %.2f\n", timeBar, wedgeHit.time, deltaTime);

            if (deltaZ <= DELTA_Z_THRESHOLD && deltaPhi <= DELTA_PHI_THRESHOLD && deltaTime <= DELTA_TIME_THRESHOLD) {
                System.out.println("  -> Qualified Cluster based on Delta thresholds.\n");
                deltaZList.add(deltaZ);
                deltaPhiList.add(deltaPhi);
                deltaTimeList.add(deltaTime);
                clusterSizes.add(3); // Two bar hits and one wedge hit
            } else {
                System.out.println("  -> Not a qualified cluster.\n");
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

        ZAndPhiClusterAnalyzer instance = new ZAndPhiClusterAnalyzer();
        instance.processEvents(reader);
        reader.close();
    }

    private void processEvents(HipoReader reader) {
        Bank atofAdcBank = new Bank(reader.getSchemaFactory().getSchema("ATOF::adc"));
        Event event = new Event();
        int eventIndex = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(atofAdcBank);

            Cluster currentCluster = new Cluster();
            int numRows = atofAdcBank.getRows();

            System.out.println("\nProcessing event " + eventIndex + " with " + numRows + " hits...");

            for (int hitIndex = 0; hitIndex < numRows; hitIndex++) {
                int sector = atofAdcBank.getInt("sector", hitIndex);
                int layer = atofAdcBank.getInt("layer", hitIndex);
                int component = atofAdcBank.getInt("component", hitIndex);
                int order = atofAdcBank.getInt("order", hitIndex);
                int adc = atofAdcBank.getInt("ADC", hitIndex);
                short pedestal = atofAdcBank.getShort("ped", hitIndex);
                double time = atofAdcBank.getFloat("time", hitIndex);

                Double z = (layer == 0) ? null : (layer == 2 || layer == 3) ? calculateZForWedge(component) : null;
                Double phi = calculatePhiForBar(component);

                EventData hit = new EventData(z, phi, time, sector, layer, component, order, adc, pedestal);
                if (layer == 0 || layer == 1) {
                    hit.printDetails("Bar");
                } else {
                    hit.printDetails("Wedge");
                }

                currentCluster.addHit(hit);
            }

            if (currentCluster.isQualifiedCluster()) {
                currentCluster.printAndRecordCluster(eventIndex);
                eventIndices.add(eventIndex);
            }

            eventIndex++;
        }

        plotDeltaVsEventIndex(deltaZList, eventIndices, "Delta Z vs Event Index", "Delta Z");
        plotDeltaVsEventIndex(deltaPhiList, eventIndices, "Delta Phi vs Event Index", "Delta Phi");
        plotDeltaVsEventIndex(deltaTimeList, eventIndices, "Delta Time vs Event Index", "Delta Time");
        plotClusterSizeVsEventIndex(clusterSizes, eventIndices);
    }

    private static double calculateZForBar(double timeLeftPMT, double timeRightPMT) {
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

    private void plotDeltaVsEventIndex(List<Double> deltaValues, List<Integer> eventIndices, String plotTitle, String yAxisLabel) {
        XYSeries series = new XYSeries(plotTitle);
        for (int i = 0; i < deltaValues.size(); i++) {
            series.add(eventIndices.get(i), deltaValues.get(i));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(plotTitle, "Event Index", yAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, plotTitle);
    }

    private void plotClusterSizeVsEventIndex(List<Integer> clusterSizes, List<Integer> eventIndices) {
        XYSeries series = new XYSeries("Cluster Size vs Event Index");
        for (int i = 0; i < clusterSizes.size(); i++) {
            series.add(eventIndices.get(i), clusterSizes.get(i));
        }
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot("Cluster Size vs Event Index", "Event Index", "Cluster Size", dataset, PlotOrientation.VERTICAL, true, true, false);
        displayChart(chart, "Cluster Size vs Event Index");
    }

    private void displayChart(JFreeChart chart, String title) {
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }
}









/*
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

public class  ZAndPhiForBarsZPositionAndTimePlotterWithDBSCAN {

    private static final int NUM_WEDGES = 10;
    private static final int NUM_BARS = 60;
    private static final double WEDGE_SPACING = 30.0;
    private static final double VELOCITY_EFF = 200.0;

    private static final double Z_THRESHOLD = 5.0;
    private static final double TIME_THRESHOLD = 1.0;

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

            // Calculation of delta values
            if (isValidCluster()) {
                double deltaZ = Math.abs(barHitLeft.zBar - wedgeHit.zWedge);
                double deltaPhi = Math.abs(barHitLeft.phiBar - wedgeHit.phiWedge);
                double deltaTime = Math.abs(barHitLeft.timeBar - wedgeHit.timeWedge);
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
        List<Double> deltaTimeList = new ArrayList<>();
        List<Double> deltaPhiList = new ArrayList<>();
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

                // Correct layer identification for wedge hits
                if (layer >= 16) {
                    layer = layer == 16 ? 2 : 3;
                }

                System.out.println("Hit -> Component: " + component + ", Layer: " + layer + ", Sector: " + sector + ", Order: " + order);

                if (component >= 1 && component <= NUM_BARS && (layer == 0 || layer == 1)) {
                    if (currentCluster.barHitLeft == null && order == 0) {
                        Double zBar = calculateZForBar(atofAdcBank, hitIndex);
                        double phiBar = calculatePhiForBar(component - 1);
                        currentCluster.addEvent(new EventData(null, zBar, null, phiBar, null, time, sector, layer, component, order, adc, pedestal));
                    } else if (currentCluster.barHitRight == null && order == 1) {
                        Double zBar = calculateZForBar(atofAdcBank, hitIndex);
                        double phiBar = calculatePhiForBar(component - 1);
                        currentCluster.addEvent(new EventData(null, zBar, null, phiBar, null, time, sector, layer, component, order, adc, pedestal));
                    }
                } else if (component >= 1 && component <= NUM_BARS && layer >= 2) {
                    int wedgeIndex = (component - 1) % NUM_WEDGES;
                    double wedgeZ = calculateZForWedge(wedgeIndex);
                    double wedgePhi = calculatePhiForBar(component - 1);
                    currentCluster.addEvent(new EventData(wedgeZ, null, wedgePhi, null, time, null, sector, layer, component, order, adc, pedestal));
                }
            }

            if (currentCluster.isValidCluster()) {
                clusters.add(currentCluster);
                System.out.println("Clusters formed in this event:");
                currentCluster.printCluster(clusters.size());

                deltaZList.add(currentCluster.barHitLeft.zBar - currentCluster.wedgeHit.zWedge);
                deltaPhiList.add(currentCluster.barHitLeft.phiBar - currentCluster.wedgeHit.phiWedge);
                deltaTimeList.add(currentCluster.barHitLeft.timeBar - currentCluster.wedgeHit.timeWedge);
                eventIndicesForDeltas.add(eventIndex);

                clusterSizes.add(currentCluster.events.size());
                clusterIndices.add(eventIndex);
            }

            eventIndex++;
        }

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

        double timeBar = Math.min(timeLeftPMT, timeRightPMT); // Use minimum time as timeBar
        double zBar = VELOCITY_EFF * (timeRightPMT - timeLeftPMT) / 2.0; // Calculate Z position based on time difference

        return timeBar; // Return timeBar as the minimum time
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

*/




/*
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

    private static final double Z_THRESHOLD = 5.0;
    private static final double TIME_THRESHOLD = 1.0;

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
                double deltaZ = Math.abs(barHitLeft.zBar - wedgeHit.zWedge);  // Z of the bar vs Z of the wedge
                double deltaPhi = Math.abs(barHitLeft.phiBar - wedgeHit.phiWedge);  // Phi of the bar vs Phi of the wedge
                double deltaTime = Math.abs(barHitLeft.timeBar - wedgeHit.timeWedge);  // Time of the bar vs Time of the wedge
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
        List<Double> deltaTimeList = new ArrayList<>();
        List<Double> deltaPhiList = new ArrayList<>();
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

                // Correct layer identification for wedge hits (layers should be 2 or 3)
                if (layer >= 16) {
                    layer = layer == 16 ? 2 : 3;
                }

                System.out.println("Hit -> Component: " + component + ", Layer: " + layer + ", Sector: " + sector + ", Order: " + order);

                // Bar hit identification: Component range is 1 to 60, and layer 0 or 1 is for bars
                if (component >= 1 && component <= NUM_BARS && (layer == 0 || layer == 1)) {
                    if (currentCluster.barHitLeft == null && order == 0) {
                        Double zBar = calculateZForBar(atofAdcBank, hitIndex);
                        double phiBar = calculatePhiForBar(component - 1);
                        currentCluster.addEvent(new EventData(null, zBar, null, phiBar, null, time, sector, layer, component, order, adc, pedestal));
                    } else if (currentCluster.barHitRight == null && order == 1) {
                        Double zBar = calculateZForBar(atofAdcBank, hitIndex);
                        double phiBar = calculatePhiForBar(component - 1);
                        currentCluster.addEvent(new EventData(null, zBar, null, phiBar, null, time, sector, layer, component, order, adc, pedestal));
                    }
                } 
                // Wedge hit identification: Component range is 1 to 60 but layer 2 or 3 is for wedges
                else if (component >= 1 && component <= NUM_BARS && layer >= 2) {
                    int wedgeIndex = (component - 1) % NUM_WEDGES;
                    double wedgeZ = calculateZForWedge(wedgeIndex);
                    double wedgePhi = calculatePhiForBar(component - 1);
                    currentCluster.addEvent(new EventData(wedgeZ, null, wedgePhi, null, time, null, sector, layer, component, order, adc, pedestal));
                }
            }

            // Only form clusters if we have two bar hits and one wedge hit
            if (currentCluster.isValidCluster()) {
                clusters.add(currentCluster);
                System.out.println("Clusters formed in this event:");
                currentCluster.printCluster(clusters.size());

                // Collect delta values
                deltaZList.add(currentCluster.barHitLeft.zBar - currentCluster.wedgeHit.zWedge);
                deltaPhiList.add(currentCluster.barHitLeft.phiBar - currentCluster.wedgeHit.phiWedge);
                deltaTimeList.add(currentCluster.barHitLeft.timeBar - currentCluster.wedgeHit.timeWedge);
                eventIndicesForDeltas.add(eventIndex);

                // Add cluster size for plotting
                clusterSizes.add(currentCluster.events.size());
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
*/

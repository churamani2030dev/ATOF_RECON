package org.example;

import java.util.ArrayList;
import java.util.List;

public class Cluster {
    private List<Hit> hits = new ArrayList<>();

    public void addHit(Hit hit) {
        hits.add(hit);
    }

    public List<Hit> getHits() {
        return hits;
    }

    public List<Integer> getHitBars() {
        List<Integer> bars = new ArrayList<>();
        for (Hit hit : hits) {
            bars.add(hit.getAdc());
        }
        return bars;
    }

    public double getAverageHitTime() {
        double totalHitTime = 0;
        for (Hit hit : hits) {
            totalHitTime += hit.getHitTime();
        }
        return hits.size() > 0 ? totalHitTime / hits.size() : 0;
    }

    public double getTotalEnergyDeposited() {
        double totalEnergy = 0;
        for (Hit hit : hits) {
            totalEnergy += hit.getEnergyDeposited();
        }
        return totalEnergy;
    }

    @Override
    public String toString() {
        return "Cluster{" +
                "hits=" + hits +
                ", hitBars=" + getHitBars() +
                ", averageHitTime=" + getAverageHitTime() +
                ", totalEnergyDeposited=" + getTotalEnergyDeposited() +
                '}';
    }
}

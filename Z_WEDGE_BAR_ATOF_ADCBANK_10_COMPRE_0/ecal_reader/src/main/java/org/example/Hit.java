package org.example;

public class Hit {
    private int adc;
    private double energyDeposited;
    private double zPosition;
    private double hitTime;

    public Hit(int adc, double energyDeposited, double zPosition, double hitTime) {
        this.adc = adc;
        this.energyDeposited = energyDeposited;
        this.zPosition = zPosition;
        this.hitTime = hitTime;
    }

    public int getAdc() {
        return adc;
    }

    public double getEnergyDeposited() {
        return energyDeposited;
    }

    public double getZPosition() {
        return zPosition;
    }

    public double getHitTime() {
        return hitTime;
    }

    @Override
    public String toString() {
        return "Hit{ADC=" + adc + ", Edep=" + energyDeposited + ", Zpos=" + zPosition + ", HitTime=" + hitTime + "}";
    }
}

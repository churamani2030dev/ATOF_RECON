package org.example;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.io.HipoReader;

public class ECALDataReader {
    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the HIPO file.");
            System.exit(1);
        }

        String hipoFilePath = args[0];

        HipoReader reader = new HipoReader();
        reader.open(hipoFilePath);

        // List available schemas
        for (Schema schema : reader.getSchemaFactory().getSchemaList()) {
            System.out.println("Schema: " + schema.getName());
        }

        // Use the correct schema names based on the list obtained above
        Schema particlesSchema = reader.getSchemaFactory().getSchema("data::event");
        Schema detectorsSchema = reader.getSchemaFactory().getSchema("data::detector");
        Schema hitsSchema = reader.getSchemaFactory().getSchema("hits::detector"); // Example schema for hits
        Schema adcSchema = reader.getSchemaFactory().getSchema("adc::detector"); // Example schema for ADC values

        if (particlesSchema == null || detectorsSchema == null || hitsSchema == null || adcSchema == null) {
            System.err.println("One or more schemas not found in the HIPO file.");
            reader.close();
            System.exit(1);
        }

        Bank particles = new Bank(particlesSchema);
        Bank detectors = new Bank(detectorsSchema);
        Bank hits = new Bank(hitsSchema);
        Bank adc = new Bank(adcSchema);
        Event event = new Event();

        int eventNumber = 0;

        while (reader.hasNext()) {
            eventNumber++;

            reader.nextEvent(event);
            event.read(particles);
            event.read(detectors);
            event.read(hits);
            event.read(adc);

            int prows = particles.getRows();
            int drows = detectors.getRows();
            int hrows = hits.getRows();
            int arows = adc.getRows();

            System.out.println("------- event # " + eventNumber);

            // Print particles information
            for (int p = 0; p < prows; p++) {
                int pid = particles.getInt("pid", p);
                float px = particles.getFloat("px", p);
                float py = particles.getFloat("py", p);
                float pz = particles.getFloat("pz", p);

                double mom = Math.sqrt(px * px + py * py + pz * pz);
                System.out.println("pid = " + pid + " mom = " + mom);

                // Print detectors information associated with the particle
                for (int d = 0; d < drows; d++) {
                    int index = (int) detectors.getByte("pindex", d);
                    if (index == p) {
                        float time = detectors.getFloat("time", d);
                        float energy = detectors.getFloat("energy", d);
                        float x = detectors.getFloat("x", d);
                        float y = detectors.getFloat("y", d);
                        float z = detectors.getFloat("z", d);
                        System.out.println(String.format("\t (%8.3f,%8.3f,%8.3f) time = %8.3f, energy = %8.3f ",
                                x, y, z, time, energy));
                    }
                }
            }

            // Print hits information
            System.out.println("Hits:");
            for (int h = 0; h < hrows; h++) {
                int hitID = hits.getInt("hitID", h);
                float hitTime = hits.getFloat("time", h);
                float hitEnergy = hits.getFloat("energy", h);
                float hitX = hits.getFloat("x", h);
                float hitY = hits.getFloat("y", h);
                float hitZ = hits.getFloat("z", h);
                System.out.println(String.format("hitID = %d, time = %8.3f, energy = %8.3f, position = (%8.3f,%8.3f,%8.3f)",
                        hitID, hitTime, hitEnergy, hitX, hitY, hitZ));
            }

            // Print ADC values
            System.out.println("ADC Values:");
            for (int a = 0; a < arows; a++) {
                int adcValue = adc.getInt("ADC", a);
                float adcTime = adc.getFloat("time", a);
                System.out.println(String.format("ADC = %d, time = %8.3f", adcValue, adcTime));
            }
        }
        reader.close();
    }
}

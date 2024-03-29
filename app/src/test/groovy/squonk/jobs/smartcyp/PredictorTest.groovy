/*
 * This Spock specification was generated by the Gradle 'init' task.
 */
package squonk.jobs.smartcyp

import org.openscience.cdk.interfaces.IAtomContainer
import org.openscience.cdk.smiles.SmilesParser
import org.openscience.cdk.DefaultChemObjectBuilder
import org.openscience.cdk.io.MDLV2000Reader
import org.openscience.cdk.io.formats.MDLV3000Format
import spock.lang.Specification

class PredictorTest extends Specification {

    def "test smiles"() {
        setup:
        def predictor = new Predictor(true, true, true, true, null, null, null, "standard")
        SmilesParser sp  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol   = sp.parseSmiles("c1ccccc1");

        when:
        def result = predictor.calculateMol(mol)

        then:
        result != null
    }

    def "test sdf v2000"() {
        setup:
        def predictor = new Predictor(true, true, true, true, null, null, null, "standard")

        when:
        predictor.run("../data/dhfr_3d-10.sdf", null)

        then:
        predictor.getCalcCount() == 10
    }

    def "test ibuprofen v2000"() {
        setup:
        def predictor = new Predictor(true, true, true, true, null, null, null, "standard")

        when:
        predictor.run("../data/Ibuprofen_V2000.sdf", null)

        then:
        predictor.getCalcCount() == 1
    }

//    def "test ibuprofen v3000"() {
//        setup:
//        def predictor = new Predictor(true, true, true, true, null, null, null, "standard")
//
//        when:
//        predictor.run("../data/Ibuprofen_V3000.sdf", null)
//
//        then:
//        predictor.getCalcCount() == 1
//    }
}
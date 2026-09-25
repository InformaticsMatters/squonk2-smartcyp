package squonk.jobs.smartcyp

import spock.lang.Specification

/**
 * Pins the wire format of the score fields.
 *
 * These assertions exist because the one-decimal-place formatting was once
 * lost: the change reached the published container image but never the
 * repository, so for a long time a build from source emitted the raw float
 * ("46.702797") where the running Job emitted "46.7". Nothing caught it -
 * PredictorTest only asserts that a result is non-null and that the molecule
 * count is right, so every score could have changed and the suite would still
 * have been green.
 *
 * The strings below are what the Data Manager writes into SMARTCyp_GEN,
 * SMARTCyp_2D6 and SMARTCyp_2C9, and what downstream consumers parse. Treat a
 * failure here as a change to the Job's output, not as a test to update.
 */
class AtomScoreSetTest extends Specification {

    private static AtomScoreSet setOf(List<Score> scores) {
        def set = new AtomScoreSet()
        scores.each { set.addScore(it) }
        set
    }

    def "V1 writes rank, atom and score to one decimal place"() {

        setup:
        def set = setOf([
                new Score(12, "N", 46.702797f, 1),
                new Score(13, "N", 48.35204f, 2),
                new Score(16, "C", 60.8f, 3),
        ])

        expect:
        set.asStringV1() == "1 N.12=46.7\n2 N.13=48.4\n3 C.16=60.8"
    }

    def "V2 writes atom index and score to one decimal place"() {

        setup:
        def set = setOf([
                new Score(12, "N", 46.702797f, 1),
                new Score(13, "N", 48.35204f, 2),
                new Score(16, "C", 60.8f, 3),
        ])

        expect:
        set.asStringV2() == "12 46.7\n13 48.4\n16 60.8"
    }

    def "a null rank leaves the rank prefix off"() {

        setup:
        def set = setOf([new Score(12, "N", 46.702797f, null)])

        expect:
        set.asStringV1() == "N.12=46.7"
    }

    def "the decimal separator does not follow the default locale"() {

        setup:
        // A locale that writes 46,7 rather than 46.7. Without Locale.ROOT in
        // AtomScoreSet the Job's output would depend on the container's LANG.
        def previous = Locale.getDefault()
        Locale.setDefault(Locale.GERMANY)
        def set = setOf([new Score(12, "N", 46.702797f, 1)])

        when:
        def result = set.asStringV1()

        then:
        result == "1 N.12=46.7"

        cleanup:
        Locale.setDefault(previous)
    }

    def "scores are rounded, not truncated"() {

        setup:
        def set = setOf([
                new Score(1, "C", 1.44f, null),
                new Score(2, "C", 1.45f, null),
                new Score(3, "C", 1.46f, null),
        ])

        expect:
        set.asStringV2() == "1 1.4\n2 1.5\n3 1.5"
    }
}

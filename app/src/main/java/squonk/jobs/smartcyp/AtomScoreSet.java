package squonk.jobs.smartcyp;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;
import java.util.Locale;
import java.util.stream.Collectors;

public class AtomScoreSet {

    /**
     * Scores are written to one decimal place.
     *
     * This is not cosmetic: it is the format the published images have emitted
     * since 2024, and downstream consumers read these fields as text. It was
     * missing from this class for a while - the change reached the container
     * image but never the repository - so a build from source silently emitted
     * the raw float ("46.702797" where the running Job said "46.7"). See
     * AtomScoreSetTest, which pins the format so that cannot happen quietly
     * again.
     *
     * Locale.ROOT, not the default locale: "%.1f" under a locale with a comma
     * decimal separator produces "46,7". The images only avoided that by
     * happening to set LANG=C.UTF-8.
     */
    private static final String SCORE_FORMAT = "%.1f";

    private List<Score> scores = new ArrayList<>();

    public List<Score> getScores() {
        return scores;
    }

    public void addScore(Score score) {
        scores.add(score);
    }

    public void sortByRank() {
        scores.sort((i, j) -> i.getRank().compareTo(j.getRank()));
    }

    public void filter(Number threshold, Integer maxRank) {

        if (threshold == null && maxRank == null) {
            return;
        }

        ListIterator<Score> iter = scores.listIterator();
        while (iter.hasNext()) {
            Score score = iter.next();
            boolean scoreFilter = threshold != null && score.getScore() > threshold.floatValue();
            boolean rankFilter = maxRank != null && score.getRank() > maxRank;
            if (scoreFilter || rankFilter) {
                iter.remove();
            }
        }
    }

    /** Score written as "1 N.12=46.7"
     * Rank-space-atomsymbol-atomnumber=score
     *
     * @return The formatted score
     */
    public String asStringV1() {
        return scores.stream().map((s) ->
                (s.getRank() == null ? "" : s.getRank() + " ")
                        + s.getAtomSymbol() + "." + s.getAtomIndex()
                        + "=" + formatScore(s)
        ).collect(Collectors.joining("\n"));
    }

    /** Score written as "12 46.7"
     * atomnumber-space-score
     *
     * @return The formatted score
     */
    public String asStringV2() {
        return scores.stream().map((s) -> {
            return s.getAtomIndex() + " " + formatScore(s);
        }).collect(Collectors.joining("\n"));
    }

    private static String formatScore(Score score) {
        return String.format(Locale.ROOT, SCORE_FORMAT, score.getScore());
    }

    public String asString(String version) {
        if (version.equals("V1")) {
            return asStringV1();
        } else if (version.equals("V2")) {
            return asStringV2();
        } else {
            throw new IllegalArgumentException("Unsupported version: " + version);
        }
    }

}
package squonk.jobs.smartcyp;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;
import java.util.stream.Collectors;

public class AtomScoreSet {

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

    /** Score written as "1 N.12=46.702797"
     * Rank-space-atomsymbol-atomnumber=score
     *
     * @return The formatted score
     */
    public String asStringV1() {
        return scores.stream().map((s) -> s.toString()).collect(Collectors.joining("\n"));
    }

    /** Score written as "12 46.702797"
     * atomnumber-space-score
     *
     * @return The formatted score
     */
    public String asStringV2() {
        return scores.stream().map((s) -> {
            return s.getAtomIndex() + " " + s.getScore();
        }).collect(Collectors.joining("\n"));
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
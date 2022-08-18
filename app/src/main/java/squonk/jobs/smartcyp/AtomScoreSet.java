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

    public String asStringV1() {
        return scores.stream().map((s) -> s.toString()).collect(Collectors.joining("\n"));
    }
}
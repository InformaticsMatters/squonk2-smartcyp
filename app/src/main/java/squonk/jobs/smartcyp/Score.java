package squonk.jobs.smartcyp;

import java.io.Serializable;

public class Score implements Serializable {

    /**
     * zero based atom index in the molecule
     */
    private final int atomIndex;
    /**
     * The atom symbol e.g. C
     */
    private final String atomSymbol;
    /**
     * The score
     */
    private final Float score;
    /**
     * Optional rank for the property, 1 being the highest. There can be ties
     */
    private final Integer rank;

    public Score(
            int atomIndex,
            String atomSymbol,
            Float score,
            Integer rank) {
        this.atomIndex = atomIndex;
        this.atomSymbol = atomSymbol;
        this.score = score;
        this.rank = rank;
    }

    public int getAtomIndex() {
        return atomIndex;
    }

    public String getAtomSymbol() {
        return atomSymbol;
    }

    public Float getScore() {
        return score;
    }

    public Integer getRank() {
        return rank;
    }

    @Override
    public String toString() {
        return (rank == null ? "" : rank + " ") + atomSymbol + "." + atomIndex + "=" + score;
    }
}
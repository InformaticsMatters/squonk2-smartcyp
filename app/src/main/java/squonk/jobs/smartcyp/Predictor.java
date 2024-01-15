/*
 * Copyright (c) 2022 Informatics Matters Ltd.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Cytochrome P450 predictions using SmartCyp
 */
package squonk.jobs.smartcyp;

import org.apache.commons.cli.*;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.DeduceBondSystemTool;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import smartcyp.MoleculeKU;
import smartcyp.SMARTSnEnergiesTable;
import squonk.jobs.util.DMLogger;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class Predictor {

    private static final Logger LOG = Logger.getLogger(Predictor.class.getName());
    private static final DMLogger DMLOG = new DMLogger();

    private static final String ORIG_ATOM_INDEX = "__OrigAtomIndex__";
    public static final String PARAM_GEN = "GEN";
    public static final String PARAM_2D6 = "2D6";
    public static final String PARAM_2C9 = "2C9";
    private static final String FIELD_PREFIX = "SMARTCyp_";
    public static final String FIELD_NAME_GEN = FIELD_PREFIX + PARAM_GEN;
    public static final String FIELD_NAME_2D6 = FIELD_PREFIX + PARAM_2D6;
    public static final String FIELD_NAME_2C9 = FIELD_PREFIX + PARAM_2C9;

    private final CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance());
    private final SMARTSnEnergiesTable _SMARTSnEnergiesTable = new SMARTSnEnergiesTable();

    private final StructureDiagramGenerator diagramGenerator = new StructureDiagramGenerator();

    private final boolean calcGeneral;
    private final boolean calc2D6;
    private final boolean calc2C9;
    private final boolean nOxidationCorrection;
    private final Integer maxRank;
    private final Float threshold;

    private Integer interval;

    private String format;

    private long calcCount = 0;


    public Predictor(boolean calcGeneral, boolean calc2D6, boolean calc2C9, boolean nOxidationCorrection,
                     Integer maxRank, Float threshold, Integer interval, String format) {
        this.calcGeneral = calcGeneral;
        this.calc2D6 = calc2D6;
        this.calc2C9 = calc2C9;
        this.nOxidationCorrection = nOxidationCorrection;
        this.maxRank = maxRank;
        this.threshold = threshold;
        this.interval = interval;
        this.format = format;
    }

    public long getCalcCount() {
        return calcCount;
    }

    public static void main(String[] args) throws Exception {

        Options options = new Options();
        options.addOption("h", "help", false, "Display help");
        options.addOption(Option.builder("i")
                .longOpt("input")
                .hasArg()
                .argName("file")
                .desc("Input file with molecules (.sdf)")
                .required()
                .build());
        options.addOption(Option.builder("o")
                .longOpt("output")
                .hasArg()
                .argName("file")
                .desc("Output file for molecules (.sdf)")
                .build());
        options.addOption(Option.builder("a")
                .longOpt("all")
                .desc("Calculate all reactivities")
                .type(Boolean.class)
                .build());
        options.addOption(Option.builder("g")
                .longOpt("calc-general")
                .desc("Calculate general reactivity")
                .type(Boolean.class)
                .build());
        options.addOption(Option.builder("d")
                .longOpt("calc-2d6")
                .desc("Calculate Cyp 2D6 reactivity")
                .type(Boolean.class)
                .build());
        options.addOption(Option.builder("c")
                .longOpt("calc-2c9")
                .desc("Calculate Cyp 2C9 reactivity")
                .type(Boolean.class)
                .build());
        options.addOption(Option.builder("e")
                .longOpt("empirical")
                .desc("Empirical N-Oxidation corrections")
                .type(Boolean.class)
                .build());
        options.addOption(Option.builder("m")
                .longOpt("max-rank")
                .hasArg()
                .argName("value")
                .desc("Maximum rank (default 3)")
                .type(Integer.class)
                .build());
        options.addOption(Option.builder("t")
                .longOpt("threshold")
                .hasArg()
                .argName("value")
                .desc("Score threshold")
                .type(Float.class)
                .build());
        options.addOption(Option.builder(null)
                .longOpt("interval")
                .hasArg()
                .argName("value")
                .desc("Reporting interval")
                .type(Integer.class)
                .build());
        options.addOption(Option.builder("f")
                .longOpt("format")
                .hasArg()
                .argName("standard/simple")
                .desc("Score format")
                .build());

        if (args.length == 0 | (args.length == 1 && ("-h".equals(args[0]) | "--help".equals(args[0])))) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(null);
            formatter.printHelp("Predictor", options);
        } else {
            CommandLineParser parser = new DefaultParser();
            CommandLine cmd = parser.parse(options, args);
            StringBuilder builder = new StringBuilder(Predictor.class.getName());
            for (String arg : args) {
                builder.append(" ").append(arg);
            }
            DMLOG.logEvent(DMLogger.Level.INFO, builder.toString());
            boolean[] methods = new boolean[3];
            if (cmd.hasOption("all")) {
                methods[0] = true;
                methods[1] = true;
                methods[2] = true;
            } else {
                methods[0] = cmd.hasOption("calc-general");
                methods[1] = cmd.hasOption("calc-2d6");
                methods[2] = cmd.hasOption("calc-2c9");
            }
            Integer maxRank = cmd.hasOption("max-rank") ? Integer.valueOf(cmd.getOptionValue("max-rank")) : null;
            Float threshold = cmd.hasOption("threshold") ? Float.valueOf(cmd.getOptionValue("threshold")) : null;
            Integer interval = cmd.hasOption("interval") ? Integer.valueOf(cmd.getOptionValue("interval")) : null;

            String format = cmd.hasOption("format") ? cmd.getOptionValue("format") : "standard";

            if (!(format.equals("simple") || format.equals("standard"))) {
                throw new IllegalArgumentException("Format option must be 'simple' or 'standard");
            }

            Predictor exec = new Predictor(methods[0], methods[1], methods[2], cmd.hasOption("empirical"),
                    maxRank, threshold, interval, format);
            exec.run(cmd.getOptionValue("input"), cmd.getOptionValue("output")
            );
        }
    }

    public void run(String input, String output) throws IOException {

        DMLOG.logEvent(DMLogger.Level.INFO, "Calculating general " + calcGeneral + ", 2D6 " + calc2D6 + ", 2C9 " + calc2C9);

        SDFWriter writer = null;
        if (output != null) {
            Path path = Paths.get(output);
            Path dir = path.getParent();
            if (dir != null) {
                Files.createDirectories(dir);
            }
            writer = new SDFWriter(new FileWriter(output));
        }
        final SDFWriter fwriter = writer;

        try {
            Stream<IAtomContainer> mols = readFile(input);
            final AtomicInteger molCount = new AtomicInteger(0);
            final AtomicInteger errorCount = new AtomicInteger(0);
            mols = mols.sequential().peek((mol) -> {
                molCount.incrementAndGet();
                if (interval != null && molCount.intValue() % interval == 0) {
                    DMLOG.logEvent(DMLogger.Level.INFO, String.format("Processed %s molecules", molCount));
                }
                IAtomContainer result = null;
                if (mol != null) {
                    try {
                        MoleculeKU molKU = calculateMol(mol);
                        writeProperties(mol, molKU);
                        if (fwriter != null) {
                            writeSDFMol(fwriter, mol);
                        }
                    } catch (CDKException | CloneNotSupportedException e) {
                        LOG.warning(String.format("Failed to process molecule %s, %s", molCount.toString(), e.getLocalizedMessage()));
                    }
                }
            }).filter((mol) -> {
                if (mol == null) {
                    errorCount.incrementAndGet();
                    return false;
                } else {
                    return true;
                }
            });

            long count = mols.count();
            calcCount += count;
            DMLOG.logEvent(DMLogger.Level.INFO, String.format("Processed %s mols, %s errors", count, errorCount));

        } finally {
            if (fwriter != null) {
                fwriter.close();
            }
        }
    }

    private void writeSDFMol(SDFWriter writer, IAtomContainer mol) throws CDKException {
        if (!GeometryTools.has2DCoordinates(mol) && !GeometryTools.has3DCoordinates(mol)) {
            if (mol instanceof IMolecule) {
                diagramGenerator.setMolecule((IMolecule)mol);
                diagramGenerator.generateCoordinates();
                mol = diagramGenerator.getMolecule();
            } else {
                LOG.warning("Molecule not an IMolecule, can't generate coordinates");
            }
        }
        mol.removeProperty("cdk:Title");
        mol.removeProperty("cdk:Remark");
        writer.write(mol);
    }

    public MoleculeKU calculateMol(IAtomContainer mol)
            throws CDKException, CloneNotSupportedException {

        IAtomContainer mol2 = prepareMolecule(mol);
        MoleculeKU moleculeKU = readMoleculeKU(mol2);
        if (moleculeKU != null) {
            moleculeKU.assignAtomEnergies(_SMARTSnEnergiesTable.getSMARTSnEnergiesTable());

            LOG.finer("\n ************** Calculating shortest distance to protonated amine **************");
            moleculeKU.calculateDist2ProtAmine();

            LOG.finer("\n ************** Calculating shortest distance to carboxylic acid **************");
            moleculeKU.calculateDist2CarboxylicAcid();

            LOG.finer("\n ************** Calculating Span2End**************");
            moleculeKU.calculateSpan2End();

            if (nOxidationCorrection) {
                LOG.finer("\n ************** Add Empirical Nitrogen Oxidation Corrections **************");
                moleculeKU.unlikelyNoxidationCorrection();
            }

            LOG.finer("\n ************** Calculating Accessabilities and Atom Scores **************");
            moleculeKU.calculateAtomAccessabilities();
            //compute 2DSASA
            double[] SASA = moleculeKU.calculateSASA();

            if (calcGeneral) moleculeKU.calculateAtomScores();
            if (calc2D6) moleculeKU.calculate2D6AtomScores();
            if (calc2C9) moleculeKU.calculate2C9AtomScores();

            LOG.finer("\n ************** Identifying, sorting and ranking C, N, P and S atoms **************");
            if (calcGeneral) {
                moleculeKU.sortAtoms();
                moleculeKU.rankAtoms();
            }
            if (calc2D6) {
                moleculeKU.sortAtoms2D6();
                moleculeKU.rankAtoms2D6();
            }
            if (calc2C9) {
                moleculeKU.sortAtoms2C9();
                moleculeKU.rankAtoms2C9();
            }
        }
        return moleculeKU;
    }

    /**
     * The input molecule must be a single fragment (not salts etc) and not have implicit hydrogens.
     *
     * @param mol
     * @return
     * @throws Exception
     */
    private MoleculeKU readMoleculeKU(IAtomContainer mol) throws CloneNotSupportedException, CDKException {

        MoleculeKU moleculeKU = null;

        if (mol != null) {

            moleculeKU = new MoleculeKU(mol, _SMARTSnEnergiesTable.getSMARTSnEnergiesTable());
//            //set the molecule title in the moleculeKU object
//            if (iAtomContainer.getProperty("SMIdbNAME") != "" && iAtomContainer.getProperty("SMIdbNAME") != null) {
//                iAtomContainer.setProperty(CDKConstants.TITLE, iAtomContainer.getProperty("SMIdbNAME"));
//            }
//            moleculeKU.setProperty(CDKConstants.TITLE, iAtomContainer.getProperty(CDKConstants.TITLE));
            moleculeKU.setProperties(mol.getProperties());
        }
        return moleculeKU;
    }

    public IAtomContainer prepareMolecule(IAtomContainer mol) throws CDKException {
        for (int atomIndex = 0; atomIndex < mol.getAtomCount(); ++atomIndex) {
            IAtom atom = mol.getAtom(atomIndex);
            atom.setProperty(ORIG_ATOM_INDEX, atomIndex);
        }
        mol = findBestFragment(mol);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
        mol = deduceBonds(mol);
        hydrogenAdder.addImplicitHydrogens(mol);
        CDKHueckelAromaticityDetector.detectAromaticity(mol);
        return mol;
    }

    public IAtomContainer deduceBonds(IAtomContainer mol) throws CDKException {
        DeduceBondSystemTool dbst = new DeduceBondSystemTool();
        mol = dbst.fixAromaticBondOrders((IMolecule) mol);
        return mol;
    }

    public IAtomContainer findBestFragment(IAtomContainer mol) {
        // Remove salts or solvents... Keep only the largest molecule
        if (!ConnectivityChecker.isConnected(mol)) {
            //System.out.println(atomContainerNr);
            IMoleculeSet fragments = ConnectivityChecker.partitionIntoMolecules(mol);

            int maxID = 0;
            int maxVal = -1;
            for (int i = 0; i < fragments.getMoleculeCount(); i++) {
                if (fragments.getMolecule(i).getAtomCount() > maxVal) {
                    maxID = i;
                    maxVal = fragments.getMolecule(i).getAtomCount();
                }
            }
            mol = fragments.getMolecule(maxID);
        }
        //end of salt removal
        return mol;
    }

    private Stream<IAtomContainer> readFile(String file) throws FileNotFoundException {
        if (file.endsWith(".sdf")) {
            return readSDF(file);
        } else if (file.endsWith(".smi")) {
            return readSmiles(file);
        }

        throw new IllegalArgumentException("Unsupported file format: " + file);
    }

    public Stream<IAtomContainer> readSDF(String file) throws FileNotFoundException {
        File sdfFile = new File(file);
        IteratingMDLReader reader = new IteratingMDLReader(
                new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance()
        );
        return createStreamFromIteratingReader(reader);
    }

    public Stream<IAtomContainer> readSmiles(String file) throws FileNotFoundException {
        File sdfFile = new File(file);
        IteratingSMILESReader reader = new IteratingSMILESReader(
                new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance()
        );
        return createStreamFromIteratingReader(reader);
    }

    private Stream<IAtomContainer> createStreamFromIteratingReader(IIteratingChemObjectReader<IAtomContainer> reader) {
        Stream<IAtomContainer> stream = StreamSupport.stream(
                        Spliterators.spliteratorUnknownSize(reader, Spliterator.ORDERED), true)
                .onClose(() -> {
                            try {
                                reader.close();
                            } catch (IOException e) {
                                LOG.log(Level.WARNING, "Failed to close reader", e);
                            }
                        }
                );
        return stream;
    }

    private void writeProperties(IAtomContainer mol, MoleculeKU moleculeKU) {

        AtomScoreSet scoresGeneral = new AtomScoreSet();
        AtomScoreSet scores2D6 = new AtomScoreSet();
        AtomScoreSet scores2C9 = new AtomScoreSet();

        for (int atomIndex = 0; atomIndex < moleculeKU.getAtomCount(); ++atomIndex) {
            IAtom currentAtom = moleculeKU.getAtom(atomIndex);
            String currentAtomType = currentAtom.getSymbol();
            if (currentAtomType.equals("C") || currentAtomType.equals("N") || currentAtomType.equals("P") || currentAtomType.equals("S")) {
                if (calcGeneral) {
                    Number score = MoleculeKU.SMARTCYP_PROPERTY.Score.get(currentAtom);
                    int rank = MoleculeKU.SMARTCYP_PROPERTY.Ranking.get(currentAtom).intValue();
                    Integer origAtomIndex = (Integer) currentAtom.getProperty(ORIG_ATOM_INDEX);
                    float f = score.floatValue();
                    scoresGeneral.addScore(new Score(origAtomIndex, currentAtomType, f, rank));
                }
                if (calc2D6) {
                    Number score = MoleculeKU.SMARTCYP_PROPERTY.Score2D6.get(currentAtom);
                    int rank = MoleculeKU.SMARTCYP_PROPERTY.Ranking2D6.get(currentAtom).intValue();
                    Integer origAtomIndex = (Integer) currentAtom.getProperty(ORIG_ATOM_INDEX);
                    float f = score.floatValue();
                    scores2D6.addScore(new Score(origAtomIndex, currentAtomType, f, rank));
                }
                if (calc2C9) {
                    Number score = MoleculeKU.SMARTCYP_PROPERTY.Score2C9.get(currentAtom);
                    int rank = MoleculeKU.SMARTCYP_PROPERTY.Ranking2C9.get(currentAtom).intValue();
                    Integer origAtomIndex = (Integer) currentAtom.getProperty(ORIG_ATOM_INDEX);
                    float f = score.floatValue();
                    scores2C9.addScore(new Score(origAtomIndex, currentAtomType, f, rank));
                }
            }
        }

        String flag = "V1";
        if (this.format.equals("simple")) {
            flag = "V2";
        }

        if (scoresGeneral.getScores().size() > 0) {
            scoresGeneral.filter(threshold, maxRank);
            scoresGeneral.sortByRank();
            mol.setProperty(FIELD_NAME_GEN, scoresGeneral.asString(flag));
        }
        if (scores2D6.getScores().size() > 0) {
            scores2D6.filter(threshold, maxRank);
            scores2D6.sortByRank();
            mol.setProperty(FIELD_NAME_2D6, scores2D6.asString(flag));
        }
        if (scores2C9.getScores().size() > 0) {
            scores2C9.filter(threshold, maxRank);
            scores2C9.sortByRank();
            mol.setProperty(FIELD_NAME_2C9, scores2C9.asString(flag));
        }
    }

}

import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { FileBlob, SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const repo = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const inputPath = path.join(repo, "results", "Supplementary_Tables_S1-S11.xlsx");
const sourcePath = process.env.ONFH_SUPPLEMENT_SOURCE ?? inputPath;
const beforeDir = path.join(repo, "qa", process.env.ONFH_RENDER_ONLY === "1" ? "supplement_final" : "supplement_before");
const renderStart = Number.parseInt(process.env.ONFH_RENDER_START ?? "0", 10);
const renderCount = Number.parseInt(process.env.ONFH_RENDER_COUNT ?? "0", 10);
const focusOnly = process.env.ONFH_FOCUS === "1";
const buildMode = process.env.ONFH_BUILD === "1";
const receptorAuditOnly = process.env.ONFH_RECEPTOR_AUDIT_ONLY === "1";
const renderOnly = process.env.ONFH_RENDER_ONLY === "1";
const renderNames = (process.env.ONFH_RENDER_SHEETS ?? "").split(",").map(x => x.trim()).filter(Boolean);

await fs.mkdir(beforeDir, { recursive: true });
const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(sourcePath));
if (buildMode) {
  const toColumnName = (n) => {
    let out = "";
    let x = n;
    while (x > 0) {
      x -= 1;
      out = String.fromCharCode(65 + (x % 26)) + out;
      x = Math.floor(x / 26);
    }
    return out;
  };

  const csvMatrix = async (relativePath) => {
    const csvText = await fs.readFile(path.join(repo, relativePath), "utf8");
    const csvWorkbook = await Workbook.fromCSV(csvText, { sheetName: "Data" });
    return csvWorkbook.worksheets.getItem("Data").getUsedRange().values;
  };

  const sectionFormat = {
    fill: "#DDEBE6",
    font: { bold: true, color: "#173B31", size: 10, name: "Carlito" },
    wrapText: true,
  };
  const headerFormat = {
    fill: "#6B8EAD",
    font: { bold: true, color: "#FFFFFF", size: 9, name: "Carlito" },
    borders: { bottom: { style: "medium", color: "#46647D" } },
    wrapText: true,
  };
  const bodyFormat = {
    font: { size: 9, name: "Carlito", color: "#1F2937" },
    borders: { bottom: { style: "thin", color: "#E2E8F0" } },
    wrapText: true,
    verticalAlignment: "top",
  };

  const writeBlock = (sheet, startRow, title, matrix) => {
    const width = matrix[0].length;
    const endColumn = toColumnName(width);
    const titleRange = sheet.getRange(`A${startRow}:${endColumn}${startRow}`);
    titleRange.format = sectionFormat;
    sheet.getRange(`A${startRow}`).values = [[title]];
    const headerRow = startRow + 1;
    sheet.getRange(`A${headerRow}:${endColumn}${headerRow}`).values = [matrix[0]];
    sheet.getRange(`A${headerRow}:${endColumn}${headerRow}`).format = headerFormat;
    if (matrix.length > 1) {
      const dataStart = headerRow + 1;
      const dataEnd = dataStart + matrix.length - 2;
      sheet.getRange(`A${dataStart}:${endColumn}${dataEnd}`).values = matrix.slice(1);
      sheet.getRange(`A${dataStart}:${endColumn}${dataEnd}`).format = bodyFormat;
    }
    return startRow + matrix.length + 3;
  };

  const lastNonEmptyRow = (sheet) => {
    const values = sheet.getUsedRange().values;
    for (let r = values.length - 1; r >= 0; r -= 1) {
      if ((values[r] ?? []).some((value) => value !== null && value !== undefined && String(value).trim() !== "")) {
        return r + 1;
      }
    }
    return 0;
  };

  const s1 = workbook.worksheets.getItem("Table S1");
  s1.getRange("A1").values = [["Table S1. Prespecified gene-set membership, overlap and selective-autophagy role/source audit"]];
  s1.getRange("A2").values = [["S1a reports frozen gene-set membership. S1b separates mechanism literature, analysis inclusion and post-analysis candidate tier; panel membership is not treated as ONFH evidence."]];
  const auditMatrix = await csvMatrix("results/selective_autophagy_receptor_audit.csv");
  const existingS1bIndex = s1.getUsedRange().values.findIndex(
    row => String(row?.[0] ?? "").startsWith("S1b. Selective-autophagy"),
  );
  const next1 = existingS1bIndex >= 0 ? existingS1bIndex + 1 : lastNonEmptyRow(s1) + 3;
  if (existingS1bIndex < 0) {
    writeBlock(s1, next1, "S1b. Selective-autophagy and mitophagy-associated gene role/source audit", auditMatrix);
  }
  const auditBodyStart = next1 + 2;
  const auditBodyEnd = auditBodyStart + auditMatrix.length - 2;
  s1.getRange(`A${next1}:I${next1}`).format.rowHeight = 60;
  s1.getRange(`A${next1 + 1}:I${next1 + 1}`).format.rowHeight = 42;
  s1.getRange(`A${auditBodyStart}:I${auditBodyEnd}`).format.rowHeight = 48;
  for (let row = auditBodyStart; row <= auditBodyEnd; row += 1) {
    if ((row - auditBodyStart) % 2 === 1) {
      s1.getRange(`A${row}:I${row}`).format.fill = "#EAF5FA";
    }
  }
  s1.getRange("A1:A220").format.columnWidth = 18;
  s1.getRange("B1:D220").format.columnWidth = 27;
  s1.getRange("E1:E220").format.columnWidth = 36;
  s1.getRange("F1:G220").format.columnWidth = 20;
  s1.getRange("H1:H220").format.columnWidth = 25;
  s1.getRange("I1:I220").format.columnWidth = 58;

  const auditSheets = [s1];
  if (!receptorAuditOnly) {
    const s7 = workbook.worksheets.getItem("Table S7");
    s7.getRange("A1").values = [["Table S7. Hallmark pathway stability, signed DoRothEA and exploratory YAP/TAZ activity audits"]];
    s7.getRange("A2").values = [["Includes leading-edge and leave-one-participant-out pathway audits, sign filtering, sampling-unit activities, tie-aware exact tests and the exploratory YAP/TAZ signature."]];
    let next7 = lastNonEmptyRow(s7) + 3;
    next7 = writeBlock(s7, next7, "S7f. Full-fit Hallmark leading-edge drivers (top 20 per selected pathway)",
      await csvMatrix("results/participant_fgsea_stability/fgsea_leading_edge_top20.csv"));
    next7 = writeBlock(s7, next7, "S7g. Leave-one-participant-out Hallmark sensitivity summary",
      await csvMatrix("results/participant_fgsea_stability/fgsea_lopo_summary.csv"));
    next7 = writeBlock(s7, next7, "S7h. Selected Hallmark membership and leading-edge overlap",
      await csvMatrix("results/participant_fgsea_stability/fgsea_selected_pathway_overlap.csv"));
    writeBlock(s7, next7, "S7i. Leave-one-participant-out method audit",
      await csvMatrix("results/participant_fgsea_stability/fgsea_lopo_method_audit.csv"));
    s7.getRange("A1:A400").format.columnWidth = 34;
    s7.getRange("B1:C400").format.columnWidth = 24;

    const s9 = workbook.worksheets.getItem("Table S9");
    s9.getRange("A1").values = [["Table S9. Serum-classifier performance, predictions, feature stability, matched permutations and paired model comparison"]];
    s9.getRange("A2").values = [["All reported predictions are out of fold. Candidate-space and Ma-comparator predictions use identical outer splits; paired uncertainty is conditional on fixed aggregated predictions and excludes model refitting."]];
    let next9 = lastNonEmptyRow(s9) + 3;
    next9 = writeBlock(s9, next9, "S9g. Aggregated paired out-of-fold predictions",
      await csvMatrix("results/figure_inputs/diag_oof_predictions_aggregated_v8.csv"));
    next9 = writeBlock(s9, next9, "S9h. Paired stratified bootstrap comparison",
      await csvMatrix("results/figure_inputs/diag_paired_model_comparison_v8.csv"));
    writeBlock(s9, next9, "S9i. Paired DeLong comparison",
      await csvMatrix("results/figure_inputs/diag_delong_model_comparison_v8.csv"));
    s9.getRange("A1:A1500").format.columnWidth = 42;
    s9.getRange("C1:D1500").format.columnWidth = 22;
    s9.getRange("G1:G1500").format.columnWidth = 50;
    s9.getRange("K1:K1500").format.columnWidth = 50;

    const s10 = workbook.worksheets.getItem("Table S10");
    let next10 = lastNonEmptyRow(s10) + 3;
    writeBlock(s10, next10, "S10i. Main-text cross-donor common-nuclear audit after mtDNA-feature exclusion",
      await csvMatrix("results/official_r_vko_no_mt_cross_donor_audit.csv"));
    s10.getRange("A1:A1000").format.columnWidth = 25;
    s10.getRange("M1:O1000").format.columnWidth = 24;
    auditSheets.push(s7, s9, s10);
  }

  for (const sheet of auditSheets) {
    sheet.showGridLines = false;
    sheet.freezePanes.freezeRows(5);
  }

  const formulaErrors = [];
  for (const sheet of auditSheets) {
    const values = sheet.getUsedRange().values;
    for (let r = 0; r < values.length; r += 1) {
      for (let c = 0; c < values[r].length; c += 1) {
        const value = values[r][c];
        if (typeof value === "string" && /^#(?:REF!|DIV\/0!|VALUE!|NAME\?|N\/A)/.test(value)) {
          formulaErrors.push(`${sheet.name}!R${r + 1}C${c + 1}:${value}`);
        }
      }
    }
  }

  const exported = await SpreadsheetFile.exportXlsx(workbook);
  await exported.save(inputPath);
  process.stdout.write(`Exported ${inputPath}\nFormula error scan: ${formulaErrors.length}\n`);
  if (formulaErrors.length) process.stdout.write(`${formulaErrors.join("\n")}\n`);
  process.exit(formulaErrors.length ? 2 : 0);
}
if (focusOnly) {
  for (const target of ["Table S7", "Table S9", "Table S10"]) {
    const sheet = workbook.worksheets.getItem(target);
    const used = sheet.getUsedRange();
    const values = used.values;
    const sectionRows = [];
    for (let i = 0; i < values.length; i += 1) {
      const first = values[i]?.[0];
      if (first !== null && first !== undefined && String(first).trim() !== "") {
        if (i < 12 || /^S\d+/i.test(String(first))) sectionRows.push([i + 1, first]);
      }
    }
    process.stdout.write(`--- ${target} section rows ---\n${JSON.stringify(sectionRows, null, 2)}\n`);
  }
  process.exit(0);
}
const sheetInfo = await workbook.inspect({
  kind: "sheet",
  include: "id,name",
  maxChars: 12000,
});
process.stdout.write(`${sheetInfo.ndjson}\n`);

const sheetsToRender = renderNames.length
  ? renderNames.map(name => workbook.worksheets.getItem(name))
  : workbook.worksheets.items.slice(renderStart, renderStart + renderCount);
const renderRanges = {
  "Table S1": "A154:I191",
  "Table S2": "A1:T33",
  "Table S3": "A1:O31",
  "Table S4": "A1:Q101",
  "Table S5": "A1:J14",
  "Table S6": "A1:E20",
  "Table S7": "A165:M279",
  "Table S8": "A1:AM40",
  "Table S9": "A1344:K1407",
  "Table S10": "A940:T968",
  "Table S11a": "A1:H36",
  "Table S11b": "A1:F74",
  "Table S11c": "A1:B80",
  "Table S11d": "A1:N60",
};
for (const sheet of sheetsToRender) {
  const safeName = sheet.name.replace(/[^A-Za-z0-9_.-]+/g, "_");
  const preview = await workbook.render({
    sheetName: sheet.name,
    range: renderRanges[sheet.name],
    autoCrop: "all",
    scale: 0.8,
    format: "png",
  });
  await fs.writeFile(
    path.join(beforeDir, `${safeName}.png`),
    new Uint8Array(await preview.arrayBuffer()),
  );
}
if (renderOnly) process.exit(0);

const overview = await workbook.inspect({
  kind: "workbook,sheet,table",
  maxChars: 20000,
  tableMaxRows: 4,
  tableMaxCols: 8,
  tableMaxCellChars: 80,
});
process.stdout.write(`${overview.ndjson}\n`);

for (const target of ["Table S7", "Table S9", "Table S10"]) {
  const region = await workbook.inspect({
    kind: "region",
    sheetId: target,
    range: "A1:T24",
    maxChars: 10000,
  });
  process.stdout.write(`--- ${target} region ---\n${region.ndjson}\n`);
  const style = await workbook.inspect({
    kind: "computedStyle",
    sheetId: target,
    range: "A1:T8",
    maxChars: 5000,
  });
  process.stdout.write(`--- ${target} style ---\n${style.ndjson}\n`);
}

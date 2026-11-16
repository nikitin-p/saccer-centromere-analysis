# Data shaping, theme and figure helpers for report.qmd.
# Inputs are the pipeline outputs (aligned.tsv, preprocess_qc.tsv) plus the
# S288C reference; every number in the report is derived here.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(ggseqlogo)
  library(ggdist)
  library(patchwork)
  library(plotly)
  library(scales)
})

CEN_LEVELS <- paste0("CEN", 1:16)
CDEI_LEN <- 8
CDEIII_LEN <- 26
BASES <- c("A", "C", "G", "T")
ELEMENTS <- c("CDEI", "CDEII", "CDEIII")

# Chart ink and palette (dataviz reference instance, light mode)
SURFACE <- "#fcfcfb"
INK <- "#0b0b0b"
INK2 <- "#52514e"
MUTED <- "#898781"
GRID <- "#e1e0d9"
AXIS <- "#c3c2b7"
BAND <- "#f0efec"
SERIES1 <- "#2a78d6"
VARIANT_COLOURS <- c(Substitution = "#2a78d6", Insertion = "#eb6834", Deletion = "#1baf7a")
DIV_LOW <- "#2a78d6"
DIV_MID <- "#f0efec"
DIV_HIGH <- "#e34948"
SEQ_LOW <- "#cde2fb"
SEQ_HIGH <- "#0d366b"

theme_cen <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "") +
    theme(
      text = element_text(colour = INK2),
      plot.title = element_text(colour = INK, face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(colour = INK2),
      axis.text = element_text(colour = MUTED),
      axis.title = element_text(colour = INK2),
      axis.line.x = element_line(colour = AXIS, linewidth = 0.3),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.text = element_text(colour = INK, face = "bold", hjust = 0),
      plot.background = element_rect(fill = SURFACE, colour = NA),
      panel.background = element_rect(fill = SURFACE, colour = NA),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(colour = INK2),
      legend.text = element_text(colour = INK2)
    )
}

save_fig <- function(p, name, width, height, dir = "figures") {
  dir.create(dir, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = SURFACE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = width, height = height, bg = SURFACE,
         device = cairo_pdf)
  invisible(p)
}

# The height is set on the widget container (not the plotly layout), and knitr
# is told not to rescale it, so the caption always sits below the plot
as_plotly <- function(p, height) {
  w <- ggplotly(p, tooltip = "text") |>
    config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d")) |>
    layout(hoverlabel = list(bgcolor = "white", font = list(color = INK)))
  w$width <- NULL
  w$height <- height
  w$sizingPolicy$knitr$figure <- FALSE
  w$sizingPolicy$defaultWidth <- "100%"
  w
}

# ---- Inputs ---------------------------------------------------------------

read_fasta <- function(path) {
  lines <- readLines(path)
  hdr <- grepl("^>", lines)
  id <- sub("^>(\\S+).*", "\\1", lines[hdr])
  seqs <- tapply(lines[!hdr], cumsum(hdr)[!hdr], paste, collapse = "")
  setNames(toupper(as.character(seqs)), id)
}

load_inputs <- function(data_dir, repo) {
  reference <- read_fasta(file.path(repo, "data/sacer_centromere.fasta"))

  aligned <- read_tsv(file.path(data_dir, "aligned.tsv"), show_col_types = FALSE) |>
    mutate(
      centromere = factor(centromere, CEN_LEVELS),
      ref_length = nchar(reference[as.character(centromere)]),
      delta_length = sequence_length - ref_length,
      blast_length = nchar(blast_sequence)
    )

  meta <- read_excel(file.path(repo, "data/41586_2018_30_MOESM3_ESM.xls"), skip = 3) |>
    select(strain = `Standardized name`, isolate = `Isolate name`,
           ecology = `Ecological origins`, geography = `Geographical origins`,
           clade = Clades, ploidy = Ploidy, zygosity = Zygosity)

  strains <- aligned |>
    distinct(sample, assembly, strain) |>
    left_join(meta, by = "strain") |>
    mutate(
      clade_no = suppressWarnings(as.integer(sub("\\..*", "", clade))),
      clade = coalesce(sub("^\\d+\\.\\s*", "", clade), "unassigned")
    ) |>
    arrange(clade_no, assembly)

  # Assemblies ordered by clade so related strains sit together in every plot
  assembly_levels <- rev(strains$assembly)
  aligned <- aligned |>
    left_join(select(strains, sample, clade, ecology), by = "sample") |>
    mutate(assembly = factor(assembly, assembly_levels))

  list(
    reference = reference,
    aligned = aligned,
    strains = strains,
    qc = read_tsv(file.path(data_dir, "preprocess_qc.tsv"), show_col_types = FALSE)
  )
}

gc_content <- function(x) {
  vapply(strsplit(x, ""), function(b) mean(b %in% c("G", "C")), numeric(1))
}

# ---- Alignment columns -> reference coordinates ---------------------------

# One row per alignment column of a CEN, with the S288C position it sits at
# (insertion columns carry the preceding reference position) and its element.
column_map <- function(ref_aligned, ref_length) {
  chars <- toupper(strsplit(ref_aligned, "")[[1]])
  is_base <- chars != "-"
  refpos <- cumsum(is_base)
  element <- case_when(
    refpos > ref_length - CDEIII_LEN ~ "CDEIII",
    refpos <= CDEI_LEN & (is_base | refpos < CDEI_LEN) ~ "CDEI",
    TRUE ~ "CDEII"
  )
  tibble(col = seq_along(chars), ref_base = chars, is_base, refpos, element)
}

cen_columns <- function(aligned) {
  aligned |>
    distinct(centromere, ref_aligned_seq, ref_length) |>
    group_by(centromere) |>
    reframe(column_map(ref_aligned_seq, ref_length))
}

# Boxes around CDEI and CDEIII in alignment-column units
cde_boxes <- function(columns) {
  columns |>
    filter(is_base, element != "CDEII") |>
    group_by(centromere, element) |>
    summarise(xmin = min(col) - 0.5, xmax = max(col) + 0.5, .groups = "drop")
}

# ---- Variants against S288C -------------------------------------------------

call_variants <- function(aligned, columns) {
  rows <- lapply(seq_len(nrow(aligned)), function(i) {
    a <- aligned[i, ]
    cols <- columns[columns$centromere == a$centromere, ]
    s <- toupper(strsplit(a$aligned_seq, "")[[1]])
    r <- cols$ref_base
    type <- case_when(
      r == "-" & s == "-" ~ NA_character_,
      r == "-" ~ "Insertion",
      s == "-" ~ "Deletion",
      r != s & r %in% BASES & s %in% BASES ~ "Substitution",
      TRUE ~ NA_character_
    )
    # Consecutive indel columns of one type form one event; every
    # substituted column is its own event
    run <- rle(ifelse(is.na(type), "none", type))
    event <- as.character(rep(seq_along(run$lengths), run$lengths))
    event <- ifelse(!is.na(type) & type == "Substitution", paste0("s", seq_along(type)), event)
    tibble(
      sample = a$sample, assembly = a$assembly, strain = a$strain, clade = a$clade,
      centromere = a$centromere, col = cols$col, refpos = cols$refpos,
      element = cols$element, ref = r, alt = s, type = type, event = event
    ) |> filter(!is.na(type))
  })
  bind_rows(rows) |>
    mutate(type = factor(type, names(VARIANT_COLOURS)),
           element = factor(element, ELEMENTS))
}

# Events: an indel spanning k columns counts once, substitutions count per column
variant_events <- function(variants) {
  variants |>
    group_by(sample, assembly, strain, clade, centromere, type, event) |>
    summarise(element = first(element), col = first(col), refpos = first(refpos),
              size = n(), ref = paste(ref, collapse = ""), alt = paste(alt, collapse = ""),
              .groups = "drop")
}

element_rates <- function(events, aligned) {
  n_asm <- aligned |> count(centromere, name = "n_assemblies")
  lengths <- aligned |>
    distinct(centromere, ref_length) |>
    mutate(CDEI = CDEI_LEN, CDEIII = CDEIII_LEN, CDEII = ref_length - CDEI_LEN - CDEIII_LEN) |>
    pivot_longer(all_of(ELEMENTS), names_to = "element", values_to = "element_length") |>
    select(-ref_length)
  kind <- events |>
    mutate(kind = if_else(type == "Substitution", "Substitutions", "Indel events")) |>
    count(centromere, element, kind, name = "n")
  expand_grid(centromere = factor(CEN_LEVELS, CEN_LEVELS), element = ELEMENTS,
              kind = c("Substitutions", "Indel events")) |>
    filter(centromere %in% aligned$centromere) |>
    left_join(mutate(kind, element = as.character(element)), by = c("centromere", "element", "kind")) |>
    left_join(lengths, by = c("centromere", "element")) |>
    left_join(n_asm, by = "centromere") |>
    mutate(n = coalesce(n, 0L),
           per_kb = 1000 * n / (n_assemblies * element_length),
           element = factor(element, ELEMENTS),
           kind = factor(kind, c("Substitutions", "Indel events")))
}

substitution_spectrum <- function(variants, aligned, columns) {
  subs <- variants |> filter(type == "Substitution") |> count(centromere, ref, alt, name = "n")
  # Opportunity: strain-sites where S288C has base b and the strain has any base
  opp <- lapply(seq_len(nrow(aligned)), function(i) {
    a <- aligned[i, ]
    cols <- columns[columns$centromere == a$centromere, ]
    s <- toupper(strsplit(a$aligned_seq, "")[[1]])
    keep <- cols$ref_base %in% BASES & s %in% BASES
    tibble(centromere = a$centromere, ref = cols$ref_base[keep])
  }) |> bind_rows() |> count(centromere, ref, name = "sites")
  classes <- expand_grid(ref = BASES, alt = BASES) |> filter(ref != alt) |>
    mutate(class = paste0(ref, "→", alt))
  list(
    by_cen = expand_grid(centromere = factor(CEN_LEVELS, CEN_LEVELS), classes) |>
      left_join(subs, by = c("centromere", "ref", "alt")) |>
      mutate(n = coalesce(n, 0L), class = factor(class, classes$class)),
    total = classes |>
      left_join(subs |> group_by(ref, alt) |> summarise(n = sum(n), .groups = "drop"), by = c("ref", "alt")) |>
      left_join(opp |> group_by(ref) |> summarise(sites = sum(sites)), by = "ref") |>
      mutate(n = coalesce(n, 0L), per_kb = 1000 * n / sites, class = factor(class, classes$class))
  )
}

# ---- Sequence logos ---------------------------------------------------------

logo_matrices <- function(aligned) {
  groups <- split(aligned$aligned_seq, droplevels(aligned$centromere))
  lapply(groups, function(seqs) {
    m <- do.call(rbind, strsplit(toupper(seqs), ""))
    pfm <- sapply(seq_len(ncol(m)), function(j) table(factor(m[, j], levels = BASES)))
    rownames(pfm) <- BASES
    pfm
  })
}

# ---- Figures in the coursework design ---------------------------------------

# Cycled hues: identity is carried by the x position, colour only separates
# neighbouring centromeres, so repeating the 8 slots is safe here
CEN_COLOURS <- setNames(rep(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                              "#e87ba4", "#008300", "#6250d6", "#e34948"), 2), CEN_LEVELS)

length_distribution <- function(aligned) {
  rng <- range(aligned$sequence_length)
  counts <- aligned |>
    group_by(centromere, sequence_length) |>
    summarise(n = n(), who = paste(sort(as.character(assembly)), collapse = ", "), .groups = "drop")
  expand_grid(centromere = factor(levels(droplevels(aligned$centromere)), CEN_LEVELS),
              sequence_length = seq(rng[1] - 2, rng[2] + 2)) |>
    left_join(counts, by = c("centromere", "sequence_length")) |>
    mutate(n = coalesce(n, 0L),
           text = paste0("<b>", centromere, ", ", sequence_length, " bp</b><br>",
                         n, if_else(n == 1, " assembly", " assemblies"),
                         if_else(n > 0, paste0("<br>", who), "")))
}

# Half-violin outlines (right) and jittered points (left) per centromere
raincloud_data <- function(aligned, value = "GC", width = 0.38) {
  cens <- levels(droplevels(aligned$centromere))
  clouds <- lapply(seq_along(cens), function(i) {
    v <- aligned[[value]][aligned$centromere == cens[i]]
    if (length(unique(v)) < 2) return(NULL)
    dens <- density(v, cut = 3, n = 256)
    tibble(centromere = cens[i], x = c(i + 0.06, i + 0.06 + width * dens$y / max(dens$y), i + 0.06),
           y = c(min(dens$x), dens$x, max(dens$x)))
  }) |> bind_rows()
  set.seed(3)
  points <- aligned |>
    mutate(i = match(as.character(centromere), cens),
           x = i - 0.16 + runif(n(), -0.07, 0.07), y = .data[[value]])
  list(clouds = mutate(clouds, centromere = factor(centromere, CEN_LEVELS)),
       points = points, cens = cens)
}

# Pipeline overview: four stages, each a column of step cards with live counts
workflow_diagram <- function(s) {
  col_x <- c(1.75, 5.05, 8.35, 11.65)
  card_w <- 3.0
  card_h <- 0.8
  card_y <- c(4.0, 2.85, 1.7)

  cards <- tribble(
    ~stage, ~row, ~title, ~body,
    1, 1, paste(s$n_asm, "T2T assemblies"), "ScRAP, O'Donnell et al. 2023\n16 chromosomes each, .fasta.gz",
    1, 2, "BLAST 16 S288C centromeres", "blastn · word size 7 · e ≤ 1e-5\none database per assembly",
    1, 3, paste(s$raw_hits, "BLAST hits"), "candidate centromere loci",
    2, 1, "One locus, one centromere", paste0("overlapping hits → best e-value\n−", s$n_cross, " cross-centromere hits"),
    2, 2, "Merge split hits", paste0(s$n_merged, " extra piece", if (s$n_merged == 1) "" else "s", " joined to its locus"),
    2, 3, paste(s$n_loci, "loci"), paste(s$n_asm, "assemblies ×", s$n_cen, "centromeres"),
    3, 1, "Cut the locus from the genome", "hit ± 150 bp, oriented like S288C\ninsertions kept",
    3, 2, "Place CDEI and CDEIII", "PWMs from the 16 S288C elements",
    3, 3, "", "",
    4, 1, "Align each centromere", "MAFFT --auto, S288C as\nreference row",
    4, 2, "Compare with S288C", "substitutions · insertions ·\ndeletions, per element",
    4, 3, "Report", "length · GC · logos · variant map\nelement rates · spectrum"
  ) |>
    mutate(x = col_x[stage], y = card_y[row])

  stages <- tibble(x = col_x, label = c("Locate", "Clean", "Annotate", "Analyse"), n = 1:4)
  down <- cards |> filter(row < 3) |> transmute(x, y = y - card_h / 2 - 0.04, yend = y - 0.27)
  across <- tibble(x = col_x[-4] + card_w / 2 + 0.06, xend = col_x[-1] - card_w / 2 - 0.06, y = 4.75)

  # CDE schematic in the third card of stage 3, proportional to S288C lengths
  sx <- col_x[3] - card_w / 2 + 0.18
  sw <- card_w - 0.36
  seg <- c(CDEI = 8, CDEII = 83, CDEIII = 26)
  ends <- sx + sw * cumsum(seg) / sum(seg)
  schematic <- tibble(element = names(seg), xmin = c(sx, ends[-3]), xmax = ends,
                      ymin = card_y[3] - 0.05, ymax = card_y[3] + 0.17,
                      fill = c(SERIES1, BAND, SERIES1))

  ggplot() +
    geom_segment(data = across, aes(x = x, xend = xend, y = y, yend = y),
                 arrow = arrow(length = unit(0.12, "in"), type = "closed"),
                 colour = AXIS, linewidth = 0.6) +
    geom_point(data = stages, aes(x - card_w / 2 + 0.17, 4.75), size = 7.5, colour = SERIES1) +
    geom_text(data = stages, aes(x - card_w / 2 + 0.17, 4.75, label = n),
              colour = "white", fontface = "bold", size = 3.6) +
    geom_text(data = stages, aes(x - card_w / 2 + 0.42, 4.75, label = label),
              hjust = 0, fontface = "bold", size = 4.6, colour = INK) +
    geom_rect(data = cards, aes(xmin = x - card_w / 2, xmax = x + card_w / 2,
                                ymin = y - card_h / 2, ymax = y + card_h / 2),
              fill = "white", colour = AXIS, linewidth = 0.35) +
    geom_rect(data = filter(cards, row == 3, stage != 3),
              aes(xmin = x - card_w / 2, xmax = x - card_w / 2 + 0.06,
                  ymin = y - card_h / 2, ymax = y + card_h / 2), fill = SERIES1) +
    geom_text(data = cards, aes(x - card_w / 2 + 0.18, y + card_h / 2 - 0.2, label = title),
              hjust = 0, vjust = 0.5, fontface = "bold", size = 3.4, colour = INK) +
    geom_text(data = cards, aes(x - card_w / 2 + 0.18, y + card_h / 2 - 0.38, label = body),
              hjust = 0, vjust = 1, size = 2.9, colour = INK2, lineheight = 0.95) +
    geom_segment(data = down, aes(x = x, xend = x, y = y, yend = yend),
                 arrow = arrow(length = unit(0.08, "in"), type = "closed"),
                 colour = AXIS, linewidth = 0.5) +
    geom_rect(data = schematic, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill)),
              colour = SURFACE, linewidth = 0.6) +
    geom_text(data = schematic, aes((xmin + xmax) / 2, ymax + 0.13, label = element),
              size = 2.7, fontface = "bold", colour = INK) +
    annotate("text", x = sx, y = card_y[3] - 0.24, hjust = 0, size = 2.8, colour = INK2,
             label = paste0("8 bp · AT-rich spacer · 26 bp (CCGAA)\n", s$n_cde_ok, "/", s$n_loci,
                            " loci with confident CDEI + CDEIII")) +
    coord_cartesian(xlim = c(0.1, 13.3), ylim = c(1.2, 5.05), expand = FALSE) +
    theme_void() +
    theme(plot.background = element_rect(fill = SURFACE, colour = NA))
}

# ---- Analysis B: CDEII in helical steps -----------------------------------
# Inputs from scripts/06_microevolution.py. Step classes are on |delta| vs the
# population median CDEII length: within <= 2, one turn 8-12, two turns 18-22,
# off-step 3-7 / 13-17, large > 22.

HELICAL_REPEAT <- 10.2
STEP_LEVELS <- c(within = "Within ±2 bp", off_step = "Off-step (3–7, 13–17)",
                 one_turn = "One turn (8–12)", two_turns = "Two turns (18–22)", large = "Large (> 22)")
STEP_COLOURS <- setNames(c(MUTED, DIV_HIGH, SERIES1, SEQ_HIGH, "#eb6834"), STEP_LEVELS)

load_microevo <- function(data_dir) {
  p <- function(f) file.path(data_dir, f)
  files <- c("cdeii_steps.tsv", "cdeii_haplotypes.tsv", "helsen_cdeii_steps.tsv", "helsen_cdeii_events.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE)
  lab <- function(x) factor(STEP_LEVELS[x], levels = STEP_LEVELS)
  list(
    steps = rd("cdeii_steps.tsv") |> mutate(class = lab(class_median),
                                            centromere = factor(centromere, levels = rev(CEN_LEVELS))),
    hap = rd("cdeii_haplotypes.tsv") |> mutate(other_assembly_cdeii = as.character(other_assembly_cdeii)),
    h_steps = rd("helsen_cdeii_steps.tsv") |> mutate(class = lab(class_median)),
    h_events = rd("helsen_cdeii_events.tsv")
  )
}

# Shaded on-step bands (8-12, 18-22 on both sides) and dashed whole-turn lines
step_bands <- function(xmax, ymin = -Inf) {
  b <- tibble(lo = c(8, 18, -12, -22), hi = c(12, 22, -8, -18)) |> filter(abs(lo) <= xmax)
  turns <- HELICAL_REPEAT * setdiff(-3:3, 0)
  list(
    geom_rect(data = b, aes(xmin = lo - 0.5, xmax = hi + 0.5, ymax = Inf), ymin = ymin,
              inherit.aes = FALSE, fill = BAND),
    geom_vline(xintercept = turns[abs(turns) <= xmax], colour = AXIS, linetype = "dashed", linewidth = 0.3),
    geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.3)
  )
}

plot_cdeii_steps <- function(steps) {
  df <- steps |>
    mutate(text = paste0("<b>", strain, "</b> ", centromere, ": CDEII ", cdeii_length,
                         " bp (median ", median_length, ", S288C ", s288c_length, ")",
                         ifelse(is.na(events) | events == "", "", paste0("<br>", events))))
  xmax <- max(24, max(abs(df$delta_median)) + 2)
  set.seed(1)
  ggplot(df, aes(delta_median, centromere, colour = class, text = text)) +
    step_bands(xmax) +
    geom_point(position = position_jitter(width = 0, height = 0.22), size = 2.2, alpha = 0.85) +
    geom_text(data = filter(df, class_median != "within"), aes(label = strain),
              nudge_y = 0.42, size = 2.8, show.legend = FALSE) +
    scale_colour_manual(values = STEP_COLOURS, name = NULL, drop = FALSE) +
    scale_x_continuous(limits = c(-xmax, xmax), breaks = seq(-30, 30, 10)) +
    labs(x = "CDEII length − population median, bp", y = NULL) +
    theme_cen()
}

plot_helsen_steps <- function(h_steps) {
  df <- h_steps |> filter(!divergent) |>
    group_by(delta_median, class) |> summarise(n = sum(n_assemblies), alleles = n(), .groups = "drop") |>
    mutate(text = paste0("Δ ", delta_median, " bp: ", n, " assemblies, ", alleles, " alleles"))
  xmax <- max(abs(df$delta_median)) + 2
  ggplot(df, aes(delta_median, n, fill = class, text = text)) +
    step_bands(xmax, ymin = log10(0.7)) +  # ymin is in transformed (log10) units
    geom_col(width = 0.8) +
    scale_fill_manual(values = STEP_COLOURS, name = NULL, drop = FALSE) +
    scale_y_log10(labels = label_comma()) +
    scale_x_continuous(breaks = seq(-30, 30, 10)) +
    labs(x = "CDEII length − population median, bp", y = "Assemblies (log scale)") +
    theme_cen()
}

# ---- Analysis C: kinetochore genes ---------------------------------------

KT_GENES <- c("CBF1", "CEP3", "NDC10", "CTF13", "SKP1", "CSE4", "SCM3", "MIF2", "OKP1", "AME1")
KT_GROUPS <- c(CBF1 = "CDEI", CEP3 = "CBF3 (CDEIII)", NDC10 = "CBF3 (CDEIII)", CTF13 = "CBF3 (CDEIII)",
               SKP1 = "CBF3 (CDEIII)", CSE4 = "Nucleosome", SCM3 = "Nucleosome",
               MIF2 = "Inner kinetochore", OKP1 = "Inner kinetochore", AME1 = "Inner kinetochore")

load_kinetochore <- function(data_dir) {
  p <- function(f) file.path(data_dir, f)
  files <- c("kinetochore_genes.tsv", "kinetochore_positions.tsv", "kinetochore_domain_summary.tsv",
             "kinetochore_strain_matrix.tsv", "kinetochore_cen_variants.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE, na = "")
  list(
    genes = rd("kinetochore_genes.tsv") |> mutate(gene = factor(gene, levels = KT_GENES)),
    positions = rd("kinetochore_positions.tsv"),
    domains = rd("kinetochore_domain_summary.tsv"),
    matrix = rd("kinetochore_strain_matrix.tsv"),
    cen = rd("kinetochore_cen_variants.tsv")
  )
}

# Strain x gene: fill = amino-acid substitutions vs S288C, label = haplotype
plot_kinetochore_matrix <- function(genes) {
  df <- genes |> filter(strain != "S288C")
  ord <- df |> filter(gene == "CBF1") |>
    arrange(as.integer(sub("H", "", haplotype)), strain) |> pull(strain)
  df <- df |>
    mutate(strain = factor(strain, levels = rev(ord)),
           group = factor(KT_GROUPS[as.character(gene)], levels = unique(KT_GROUPS)),
           fs = status != "intact",
           label = case_when(fs ~ "fs", haplotype == "ref" ~ "·", TRUE ~ haplotype),
           text = paste0("<b>", strain, "</b> ", gene, " · ", haplotype, " (", status, ")",
                         "<br>", n_substitutions, " substitutions",
                         ifelse(is.na(variants) | variants == "", "",
                                paste0("<br>", gsub(",", ", ", variants)))))
  ggplot(df, aes(gene, strain, text = text)) +
    geom_tile(aes(fill = n_substitutions), colour = SURFACE, linewidth = 0.8) +
    geom_tile(data = filter(df, fs), fill = DIV_HIGH, colour = SURFACE, linewidth = 0.8) +
    geom_text(aes(label = label, colour = n_substitutions > 12 | fs), size = 3.2) +
    scale_fill_gradient(low = BAND, high = SEQ_HIGH, name = "Substitutions vs S288C") +
    scale_colour_manual(values = c(`FALSE` = INK, `TRUE` = "white"), guide = "none") +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL) +
    theme_cen() + theme(panel.grid = element_blank(), axis.line.x = element_blank())
}

# Substitutions per 100 aa in each annotated region, S288C-private positions excluded
plot_domain_rates <- function(domains) {
  df <- domains |>
    mutate(gene = factor(gene, levels = KT_GENES),
           annotated = if_else(domain == "unannotated", "Unannotated remainder", "Annotated region"),
           domain = factor(domain, levels = rev(c(setdiff(unique(domain), "unannotated"), "unannotated"))),
           text = paste0("<b>", gene, "</b> ", domain, ": ", n_substituted_excl_ref_private,
                         " substituted of ", length, " aa (", n_intact, " intact alleles)"))
  ggplot(df, aes(subst_per_100aa, domain, colour = annotated, text = text)) +
    geom_segment(aes(x = 0, xend = subst_per_100aa, yend = domain), linewidth = 0.4) +
    geom_point(aes(size = length)) +
    scale_colour_manual(values = c(`Annotated region` = SERIES1, `Unannotated remainder` = MUTED), name = NULL) +
    scale_size_area(max_size = 4, guide = "none") +
    facet_grid(gene ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "Variable positions per 100 aa", y = NULL) +
    theme_cen() + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
}

# Analysis G: aneuploidy vs CEN allele class (Peter 2018 748-strain BLAST set) -------------

ANEU_PREDICTORS <- c(non_modal = "Not the major allele", length_variant = "Length differs from mode")

load_aneuploidy <- function(dir) {
  p <- function(f) file.path(dir, f)
  files <- c("aneuploidy_association.tsv", "aneuploidy_counts.tsv", "cen_calls.tsv", "aneuploidy_parsed.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE, na = c("", "NA"))
  list(
    assoc = rd("aneuploidy_association.tsv") |>
      mutate(centromere = factor(centromere, levels = CEN_LEVELS)),
    counts = rd("aneuploidy_counts.tsv"),
    calls = rd("cen_calls.tsv"),
    meta = rd("aneuploidy_parsed.tsv")
  )
}

plot_aneuploidy_or <- function(assoc) {
  df <- assoc |>
    filter(tested, predictor %in% names(ANEU_PREDICTORS)) |>
    mutate(predictor = factor(ANEU_PREDICTORS[predictor], levels = ANEU_PREDICTORS),
           centromere = factor(centromere, levels = rev(CEN_LEVELS)),
           estimable = !separation & is.finite(clogit_or),
           text = paste0("<b>", centromere, "</b> ", predictor,
                         "<br>aneuploid: ", aneu_variant, " of ", n_aneuploid, " carry the variant",
                         "<br>euploid: ", eu_variant, " of ", eu_variant + eu_ref,
                         "<br>OR (within clade) ", ifelse(estimable, signif(clogit_or, 2), "n/e"),
                         "<br>permutation p ", signif(perm_p, 2), ", q ", signif(perm_q, 2)))
  est <- filter(df, estimable)
  ggplot(df, aes(y = centromere, text = text)) +
    geom_vline(xintercept = 1, colour = AXIS) +
    geom_linerange(data = est, aes(xmin = pmax(clogit_lo, 0.01), xmax = pmin(clogit_hi, 100)),
                   colour = SERIES1, linewidth = 0.5) +
    geom_point(data = est, aes(x = clogit_or), colour = SERIES1, size = 2) +
    geom_point(data = filter(df, !estimable), aes(x = 1), shape = 4, colour = MUTED, size = 2) +
    scale_x_log10(limits = c(0.01, 100), breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100"),
                  oob = squish) +
    facet_wrap(~predictor, nrow = 1) +
    labs(x = "Odds ratio for aneuploidy of the chromosome, within clade (95 % CI)", y = NULL) +
    theme_cen()
}

# Analysis D: Ty elements and solo LTRs around centromeres ---------------------------------

TE_FAMILIES <- c("Ty1", "Ty2", "Ty3", "Ty4", "Ty5")
TE_COLOURS <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"), TE_FAMILIES)
TE_CLASS <- c(full = "Full element", truncated = "Truncated element", solo_LTR = "Solo LTR")
NULL_LABELS <- c(`tRNA-matched` = "tRNA-matched null", uniform = "Uniform null")

load_transposons <- function(data_dir) {
  p <- function(f) file.path(data_dir, f)
  files <- c("te_elements.tsv", "te_near_cen.tsv", "trna_near_cen.tsv", "te_cen_summary.tsv",
             "te_null.tsv", "te_sites.tsv", "te_s288c_check.tsv", "te_loegler_s5.tsv", "te_flank_calls.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE, na = c("", "NA"))
  list(
    elements = rd("te_elements.tsv"),
    near = rd("te_near_cen.tsv") |> mutate(centromere = factor(centromere, levels = CEN_LEVELS)),
    trna = rd("trna_near_cen.tsv") |> mutate(centromere = factor(centromere, levels = CEN_LEVELS)),
    summary = rd("te_cen_summary.tsv"),
    null = rd("te_null.tsv"),
    sites = rd("te_sites.tsv"),
    check = rd("te_s288c_check.tsv"),
    s5 = read_tsv(p("te_loegler_s5.tsv"), show_col_types = FALSE, col_types = cols(region = "c")),
    calls = rd("te_flank_calls.tsv")
  )
}

# Elements within +-20 kb of each CEN, per genome; CDEI side to the left
plot_te_map <- function(near, trna, genomes) {
  lev <- rev(c(setdiff(sort(genomes), "S288C"), "S288C"))
  prep <- function(df) df |> mutate(genome = factor(genome, levels = lev),
                                    lo = rel_lo / 1000, hi = rel_hi / 1000, mid = (lo + hi) / 2)
  te <- prep(near) |>
    mutate(class = factor(TE_CLASS[class], levels = TE_CLASS),
           family = factor(family, levels = TE_FAMILIES),
           text = paste0("<b>", genome, "</b> ", centromere, ": ", family, " ", class,
                         "<br>", start, "–", end, " (", end - start + 1, " bp), ",
                         distance, " bp from the CEN, ", side, " side, ", rel_strand, " strand"))
  tr <- prep(trna) |> mutate(text = paste0("<b>", genome, "</b> ", centromere, ": tRNA, ",
                                           distance, " bp from the CEN"))
  full <- filter(te, class != "Solo LTR")
  ggplot(mapping = aes(y = genome)) +
    geom_vline(xintercept = 0, colour = AXIS) +
    geom_point(data = tr, aes(x = mid, text = text), shape = 124, size = 2.2, colour = MUTED) +
    geom_segment(data = full, aes(x = lo, xend = hi, yend = genome, colour = family, text = text),
                 linewidth = 2.2) +
    geom_point(data = filter(te, class == "Solo LTR"), aes(x = mid, colour = family, text = text),
               size = 1.8) +
    scale_colour_manual(values = TE_COLOURS, name = NULL) +
    scale_x_continuous(limits = c(-20.5, 20.5), breaks = c(-20, -10, 0, 10, 20)) +
    facet_wrap(~centromere, ncol = 4, drop = FALSE) +
    labs(x = "kb from the centromere (CDEI side negative)", y = NULL) +
    theme_cen(10) + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7))
}

# Observed elements in each window over the 16 CENs vs the null 5-95 % range, per genome
plot_te_null <- function(null) {
  df <- null |> filter(genome != "all") |>
    mutate(null = factor(NULL_LABELS[null], levels = NULL_LABELS),
           window = factor(paste0("±", window_kb, " kb"), levels = paste0("±", sort(unique(window_kb)), " kb")),
           genome = factor(genome, levels = rev(sort(unique(genome)))),
           text = paste0("<b>", genome, "</b> ", window, ", ", null,
                         "<br>observed ", observed, ", null mean ", round(null_mean, 1),
                         " (5–95 %: ", null_q05, "–", null_q95, "), p ", signif(p, 2)))
  ggplot(df, aes(y = genome, text = text)) +
    geom_linerange(aes(xmin = null_q05, xmax = null_q95, colour = null), linewidth = 2.4, alpha = 0.35,
                   position = position_dodge(width = 0.6)) +
    geom_point(aes(x = null_mean, colour = null), shape = 3, size = 1.6, position = position_dodge(width = 0.6)) +
    geom_point(aes(x = observed), colour = INK, size = 2) +
    scale_colour_manual(values = c(`tRNA-matched null` = SERIES1, `Uniform null` = MUTED), name = NULL) +
    facet_wrap(~window, nrow = 1, scales = "free_x") +
    labs(x = "Elements within the window, summed over the 16 centromeres (dot: observed; bar: null 5–95 %, +: null mean)",
         y = NULL) +
    theme_cen()
}

# ---------------------------------------------------------------- analysis E: CDE-like motifs in LTRs
LTR_SETS <- c(LTR = "Ty LTRs", shuffled = "Dinucleotide-shuffled LTRs", genomic = "Random genomic windows")
LTR_SET_COLOURS <- setNames(c(SERIES1, "#eb6834", MUTED), LTR_SETS)

load_ltr_cde <- function(data_dir) {
  p <- function(f) file.path(data_dir, f)
  files <- c("ltr_cde_scores.tsv", "ltr_cde_control.tsv", "ltr_cde_poscontrol.tsv", "ltr_cde_candidates.tsv",
             "ltr_cde_score_dist.tsv", "ltr_cde_genome_units.tsv", "ltr_cde_thresholds.tsv", "filtered.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE, na = c("", "NA"))
  cand <- if (length(readLines(p("ltr_cde_candidates.tsv"), n = 1))) rd("ltr_cde_candidates.tsv") else tibble()
  list(scores = rd("ltr_cde_scores.tsv"), control = rd("ltr_cde_control.tsv"),
       pos = rd("ltr_cde_poscontrol.tsv"), candidates = cand, dist = rd("ltr_cde_score_dist.tsv"),
       units = rd("ltr_cde_genome_units.tsv"), thr = rd("ltr_cde_thresholds.tsv"),
       cen = rd("filtered.tsv") |> select(strain, centromere, cdei_score, cdeiii_score))
}

# Best CDEI / CDEIII score per sequence: LTRs vs controls; dashed = threshold, ticks = real CENs
plot_ltr_cde_scores <- function(l) {
  long <- function(df, cols) df |> pivot_longer(all_of(cols), names_to = "motif", values_to = "score") |>
    mutate(motif = factor(c(best_cdei = "CDEI", best_cdeiii = "CDEIII",
                            cdei_score = "CDEI", cdeiii_score = "CDEIII")[motif], levels = c("CDEI", "CDEIII")))
  d <- long(l$dist, c("best_cdei", "best_cdeiii")) |> mutate(set = factor(LTR_SETS[set], levels = LTR_SETS))
  cen <- long(l$cen, c("cdei_score", "cdeiii_score")) |>
    mutate(text = paste0("<b>", strain, "</b> ", centromere, ": score ", score))
  thr <- tibble(motif = factor(c("CDEI", "CDEIII"), levels = c("CDEI", "CDEIII")),
                x = c(l$thr$cdei_min, l$thr$cdeiii_min))
  ggplot(d, aes(score, colour = set)) +
    geom_density(linewidth = 0.8, adjust = 1.2) +
    geom_vline(data = thr, aes(xintercept = x), linetype = "dashed", colour = INK2) +
    geom_rug(data = cen, aes(x = score, text = text), inherit.aes = FALSE, colour = INK, alpha = 0.5,
             length = grid::unit(0.05, "npc")) +
    scale_colour_manual(values = LTR_SET_COLOURS, name = NULL) +
    facet_wrap(~motif, scales = "free", nrow = 1) +
    labs(x = "Best PWM score in the sequence (log2-odds)", y = "Density") +
    theme_cen() + theme(legend.position = "bottom")
}

# ---- Analysis F: diversity at the CEN vs its flanks --------------------------

load_diversity <- function(data_dir) {
  p <- function(f) file.path(data_dir, f)
  files <- c("diversity_windows.tsv", "diversity_summary.tsv", "diversity_crosscheck.tsv",
             "l3_loci.tsv", "flanks/windows.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE, na = c("", "NA"))
  list(win = rd("diversity_windows.tsv"), summary = rd("diversity_summary.tsv"),
       cross = rd("diversity_crosscheck.tsv"), l3 = rd("l3_loci.tsv"), windows = rd("flanks/windows.tsv"))
}

# pi along CEN +- 2 kb: x = bp from CDEI start for the left flank, from CDEIII end for the
# right flank, with the CEN drawn to a fixed 300-bp width in between so CDEs are readable
plot_pi_profile <- function(div, cen_width = 300) {
  w <- div$win |> mutate(centromere = factor(centromere, levels = CEN_LEVELS))
  fl <- w |> filter(!region %in% ELEMENTS) |>
    mutate(x = ifelse(offset < 0, (offset + 0.5) * 100, cen_width + (offset - 0.5) * 100))
  el_x <- c(CDEI = 0.15, CDEII = 0.5, CDEIII = 0.85) * cen_width
  el <- w |> filter(region %in% ELEMENTS) |> mutate(x = el_x[region], region = factor(region, levels = ELEMENTS))
  med <- bind_rows(fl, el) |> group_by(x) |> summarise(pi = median(pi, na.rm = TRUE), .groups = "drop")
  fl <- fl |> mutate(text = paste0("<b>", centromere, "</b> ", region, ": π ", signif(pi, 2),
                                   " (subst. ", signif(pi_sub, 2), ", indel ", signif(indel_rate, 2), ")"))
  el <- el |> mutate(text = paste0("<b>", centromere, "</b> ", region, ": π ", signif(pi, 2),
                                   " (subst. ", signif(pi_sub, 2), ", indel ", signif(indel_rate, 2), ")"))
  brks <- c(-2000, -1000, 0, cen_width, cen_width + 1000, cen_width + 2000)
  ggplot() +
    annotate("rect", xmin = 0, xmax = cen_width, ymin = -Inf, ymax = Inf, fill = BAND) +
    geom_line(data = fl |> filter(offset < 0), aes(x, pi, group = centromere), colour = AXIS, linewidth = 0.3) +
    geom_line(data = fl |> filter(offset > 0), aes(x, pi, group = centromere), colour = AXIS, linewidth = 0.3) +
    geom_point(data = fl, aes(x, pi, text = text), colour = AXIS, size = 0.6) +
    geom_line(data = med |> filter(x < 0), aes(x, pi), colour = INK, linewidth = 0.9) +
    geom_line(data = med |> filter(x > cen_width), aes(x, pi), colour = INK, linewidth = 0.9) +
    geom_point(data = el, aes(x, pi, colour = region, text = text),
               position = position_jitter(width = 12, height = 0, seed = 1), size = 1.6, alpha = 0.8) +
    geom_point(data = med |> filter(x >= 0, x <= cen_width), aes(x, pi), shape = 95, size = 9, colour = INK) +
    scale_colour_manual(values = c(CDEI = "#2a78d6", CDEII = "#eb6834", CDEIII = "#1baf7a"), name = NULL) +
    scale_x_continuous(breaks = brks, labels = c("−2 kb", "−1 kb", "CDEI", "CDEIII", "+1 kb", "+2 kb")) +
    scale_y_sqrt() +
    labs(x = NULL, y = "π (pairwise differences per site, √ scale)") +
    theme_cen()
}

# ---- Analysis F2: CDEII A-tract grammar and DNA shape --------------------------

load_grammar <- function(data_dir) {
  p <- function(f) file.path(data_dir, "cdeii_grammar", f)
  files <- c("atracts.tsv", "grammar_null.tsv", "shape_profiles.tsv", "shape_null.tsv",
             "shape_variants.tsv", "shape_variants_summary.tsv", "shape_events.tsv", "shape_events_summary.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE, na = c("", "NA"))
  lapply(setNames(files, sub("\\.tsv$", "", files)), rd)
}

NULL_COLOURS <- c(shuffle = "#2a78d6", genomic = "#eb6834")
GRAMMAR_LABELS <- c(n_tracts = "A-tracts per CDEII", max_tract = "Longest A-tract (bp)",
                    frac_in_tracts = "Fraction of bp in A-tracts", phase = "Phasing score",
                    power = "AA/TT power at 10–11 bp")

# Observed statistic per CEN (black) against the null median and 5-95 % range of each null
plot_grammar_null <- function(gr, stats = names(GRAMMAR_LABELS)) {
  g <- gr$grammar_null |> filter(set == "modal", centromere != "all", statistic %in% stats) |>
    mutate(centromere = factor(centromere, levels = rev(CEN_LEVELS)),
           statistic = factor(statistic, levels = stats, labels = GRAMMAR_LABELS[stats]),
           text = paste0("<b>", centromere, "</b> ", null, ": obs ", signif(obs, 3), ", null median ",
                         signif(null_median, 3), " [", signif(null_q05, 3), "–", signif(null_q95, 3),
                         "], p ", signif(p, 2)))
  dodge <- position_dodge(width = 0.6)
  ggplot(g, aes(y = centromere)) +
    geom_linerange(aes(xmin = null_q05, xmax = null_q95, colour = null), position = dodge,
                   linewidth = 1.2, alpha = 0.5) +
    geom_point(aes(x = null_median, colour = null, text = text), position = dodge, size = 1.4) +
    geom_point(data = distinct(g, centromere, statistic, obs), aes(x = obs), colour = INK, shape = 18, size = 2.4) +
    facet_wrap(~statistic, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = NULL_COLOURS, name = "Null (median, 5–95 %)",
                        labels = c(shuffle = "dinucleotide shuffles", genomic = "AT-matched S288C windows")) +
    labs(x = NULL, y = NULL) + theme_cen() + theme(legend.position = "bottom")
}

# Shape along CDEII: one grey line per modal allele, black = median over the 16
plot_shape_profiles <- function(gr, features = c("MGW", "ProT", "Roll", "HelT")) {
  units <- c(MGW = "MGW (Å)", ProT = "ProT (°)", Roll = "Roll (°)", HelT = "HelT (°)")
  pr <- gr$shape_profiles |> filter(set == "modal", feature %in% features) |>
    mutate(feature = factor(feature, levels = features, labels = units[features]))
  med <- pr |> group_by(feature, pos) |> filter(n() >= 8) |>
    summarise(value = median(value), .groups = "drop")
  ggplot(pr, aes(pos, value)) +
    geom_line(aes(group = centromere, text = centromere), colour = AXIS, linewidth = 0.3) +
    geom_line(data = med, colour = INK, linewidth = 0.8) +
    facet_wrap(~feature, scales = "free_y", ncol = 1) +
    labs(x = "Position in CDEII (bp from CDEI)", y = NULL) + theme_cen()
}

# Percentile of natural substitutions among null mutations, per feature and test level
plot_shape_variants <- function(gr) {
  ev <- gr$shape_events |>
    pivot_longer(matches("^pct_.*_(base|trinuc)$"), names_to = c("feature", "null"),
                 names_pattern = "pct_(.*)_(base|trinuc)", values_to = "pct") |>
    filter(!is.na(pct)) |>
    mutate(level = ifelse(null == "base", "Events:\nsame base", "Events:\nsame trinuc."))
  al <- gr$shape_variants |> filter(set == "Helsen") |>
    pivot_longer(starts_with("pct_"), names_to = "feature", names_prefix = "pct_", values_to = "pct") |>
    mutate(level = "Alleles:\nspectrum null")
  d <- bind_rows(al |> select(feature, level, pct), ev |> select(feature, level, pct)) |>
    mutate(level = factor(level, levels = c("Alleles:\nspectrum null", "Events:\nsame base", "Events:\nsame trinuc.")),
           feature = factor(feature, levels = c("MGW", "ProT", "Roll", "HelT")))
  ggplot(d, aes(pct)) +
    geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.1), fill = MUTED, colour = SURFACE) +
    geom_hline(yintercept = 1, colour = INK, linetype = "dashed", linewidth = 0.4) +
    facet_grid(level ~ feature) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "Percentile of the natural shape change among null mutations (0 = smallest)",
         y = "Density (1 = no preference)") + theme_cen()
}

# ---- Analysis H: centromere tree vs genome tree ------------------------------

load_trees <- function(data_dir) {
  p <- function(f) file.path(data_dir, "trees", f)
  files <- c("tree_compare.tsv", "strain_discordance.tsv", "cen_contribution.tsv", "swap_check.tsv",
             "swap_scan.tsv", "genome_core_stats.tsv", "cen_nj.nwk", "genome_nj.nwk", "cen_blocks.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE, na = c("", "NA"))
  list(compare = rd("tree_compare.tsv"), disc = rd("strain_discordance.tsv"),
       contrib = rd("cen_contribution.tsv"), swap = rd("swap_check.tsv"), scan = rd("swap_scan.tsv"),
       core = rd("genome_core_stats.tsv"), blocks = rd("cen_blocks.tsv"),
       cen_tree = ape::read.tree(p("cen_nj.nwk")), genome_tree = ape::read.tree(p("genome_nj.nwk")))
}

# Tanglegram (base graphics): genome NJ left, CEN NJ right, matching tips joined
plot_tanglegram <- function(tr, name = NULL, dir = "figures") {
  draw <- function() {
    tips <- tr$genome_tree$tip.label
    a <- cbind(tips, tips)
    op <- par(bg = SURFACE, col = INK2, mar = c(1, 1, 2, 1))
    on.exit(par(op))
    ape::cophyloplot(ape::ladderize(tr$genome_tree), ape::ladderize(tr$cen_tree), assoc = a,
                     length.line = 2, space = 40, gap = 2, col = MUTED, lwd = 1)
    u <- par("usr")
    mtext("Genome (core SNPs)", side = 3, at = u[1], adj = 0, line = 0.5, col = INK, cex = 0.9)
    mtext("Centromeres (16 concatenated)", side = 3, at = u[2], adj = 1, line = 0.5, col = INK, cex = 0.9)
  }
  if (!is.null(name)) {
    dir.create(dir, showWarnings = FALSE)
    png(file.path(dir, paste0(name, ".png")), width = 9, height = 5, units = "in", res = 300, bg = SURFACE)
    draw(); dev.off()
    cairo_pdf(file.path(dir, paste0(name, ".pdf")), width = 9, height = 5, bg = SURFACE)
    draw(); dev.off()
  }
  draw()
}

# CEN alleles vs Loegler 2025 phenotypes (EMMAX LMM, LOCO kinship) -------------------------

PHENO_FAMILIES <- c(primary_growth = "Segregation-stress growth", primary_cis = "Cis genes (±10 kb)",
                    primary_kinetochore = "Kinetochore genes", secondary = "All growth traits")
PHENO_COLOURS <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#898781"), PHENO_FAMILIES)
PHENO_PREDICTORS <- c(non_modal = "Not the major allele", cdeii_delta = "CDEII Δ length",
                      off_step = "Off a helical step", rare_burden = "Rare-allele burden")

load_phenotypes <- function(dir) {
  p <- function(f) file.path(dir, f)
  files <- c("assoc_primary.tsv", "assoc_secondary.tsv", "positive_controls.tsv",
             "trait_samples.tsv", "cen_predictors.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE, na = c("", "NA"))
  list(primary = rd("assoc_primary.tsv"), secondary = rd("assoc_secondary.tsv"),
       controls = rd("positive_controls.tsv"), samples = rd("trait_samples.tsv"),
       predictors = rd("cen_predictors.tsv"))
}

plot_pheno_qq <- function(primary, secondary) {
  df <- bind_rows(primary, secondary) |>
    filter(!is.na(p)) |>
    group_by(family) |>
    arrange(p, .by_group = TRUE) |>
    mutate(expected = -log10(ppoints(n())), observed = -log10(p),
           family = factor(PHENO_FAMILIES[family], levels = PHENO_FAMILIES),
           text = paste0("<b>", trait, "</b> ", centromere, " · ", PHENO_PREDICTORS[predictor],
                         "<br>β ", signif(beta, 2), " (SE ", signif(se, 2), "), n ", n,
                         "<br>p ", signif(p, 2), ", q ", signif(q, 2),
                         ifelse(is.na(p_local), "", paste0("<br>p | local haplotype ", signif(p_local, 2))))) |>
    ungroup()
  lim <- max(c(df$expected, df$observed))
  ggplot(df, aes(expected, observed, colour = family, text = text)) +
    geom_abline(slope = 1, intercept = 0, colour = AXIS) +
    geom_point(size = 1.1, alpha = 0.8) +
    scale_colour_manual(values = PHENO_COLOURS, name = NULL) +
    coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
    facet_wrap(~family, nrow = 1) +
    labs(x = "Expected −log10 p", y = "Observed −log10 p") +
    theme_cen() + theme(legend.position = "none")
}

# ---------------------------------------------------------------- locators L1-L4 (locator spec phase 4)
load_locators <- function(data_dir) {
  p <- function(...) file.path(data_dir, ...)
  files <- c("loci_reconciled.tsv", "l2_loci.tsv", "l3_loci.tsv", "hmm/hmm_loci.tsv",
             "hmm_from_s288c/hmm_rounds.tsv", "cen_like_elements.tsv", "cen_like_nearest.tsv",
             "locator_synthetic/locator_synthetic_summary.tsv", "locator_synthetic/locator_synthetic.tsv")
  if (!all(file.exists(p(files)))) return(NULL)
  rd <- function(f) read_tsv(p(f), show_col_types = FALSE, na = c("", "NA"))
  list(rec = rd("loci_reconciled.tsv"), l2 = rd("l2_loci.tsv"), l3 = rd("l3_loci.tsv"),
       l4 = rd("hmm/hmm_loci.tsv"), rounds = rd("hmm_from_s288c/hmm_rounds.tsv"),
       cen_like = rd("cen_like_elements.tsv"), nearest = rd("cen_like_nearest.tsv"),
       syn = rd("locator_synthetic/locator_synthetic_summary.tsv"),
       syn_detail = rd("locator_synthetic/locator_synthetic.tsv"))
}

# Per locator: loci found, overlapping L1, identical boundaries to L1
locator_vs_l1 <- function(lc) {
  row <- function(name, found, agree, exact, n) tibble(Locator = name, Found = found, `Overlaps L1` = agree,
                                                      `Same boundaries as L1` = exact, Of = n)
  l3 <- lc$l3 |> filter(status %in% c("lifted", "bridged"))
  l4 <- lc$l4 |> filter(status %in% c("ok", "insertion", "cde", "length"))
  n <- nrow(lc$rec)
  bind_rows(
    row("L1  BLAST of the S288C CEN", sum(grepl("L1", lc$rec$support)), NA, NA, n),
    row("L2  BLAST of the 2-kb flanks", sum(lc$l2$status == "ok"), sum(lc$l2$agree_l1, na.rm = TRUE),
        sum(lc$l2$d_start_l1 == 0 & lc$l2$d_end_l1 == 0, na.rm = TRUE), n),
    row("L3  minimap2 synteny + liftover", nrow(l3), sum(l3$agree, na.rm = TRUE),
        sum(l3$status == "lifted" & l3$d_start == 0 & l3$d_end == 0, na.rm = TRUE), n),
    row("L4  per-CEN HMM (nhmmer)", nrow(l4),   # agree_l1 is NA for status cde: recompute
        sum(l4$contig == l4$l1_contig & l4$hmm_start <= l4$l1_end & l4$hmm_end >= l4$l1_start, na.rm = TRUE),
        sum(l4$d_start_l1 == 0 & l4$d_end_l1 == 0, na.rm = TRUE), n)
  ) |> mutate(across(c(`Overlaps L1`, `Same boundaries as L1`), ~ ifelse(is.na(.x), "–", as.character(.x))))
}

# Profiles grown from S288C alone: sequences per profile and loci accepted, per round
plot_hmm_growth <- function(rounds) {
  long <- rounds |>
    transmute(round, centromere, `Sequences in the profile` = n_profile,
              `Assemblies found (status ok, of 10)` = found_ok, `S288C self-score, bits` = s288c_bits) |>
    pivot_longer(-c(round, centromere), names_to = "panel") |>
    mutate(panel = factor(panel, levels = c("Sequences in the profile", "Assemblies found (status ok, of 10)",
                                            "S288C self-score, bits")))
  avg <- long |> group_by(round, panel) |> summarise(value = mean(value), .groups = "drop")
  ggplot(long, aes(round, value)) +
    geom_line(aes(group = centromere), colour = AXIS, linewidth = 0.4) +
    geom_point(colour = MUTED, size = 0.9) +
    geom_line(data = avg, colour = SERIES1, linewidth = 1.1) +
    geom_point(data = avg, colour = SERIES1, size = 2) +
    facet_wrap(~panel, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(long$round))) +
    labs(x = "Round (round 0 = the S288C sequence alone)", y = NULL) +
    theme_cen()
}

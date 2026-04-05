# OSD-157 / GLDS-157: Exploracao Sistematica da Resposta Transcriptomica a Radiacao com Foco Neurodegenerativo

Autor: Arthur Pelozo  
Data: Abril 2026  
Plataforma: Agilent Whole Human Genome Microarray (Single-Channel)

## Visao geral do projeto

Este repositorio documenta uma analise exploratoria e explicativa do efeito da radiacao ionizante (0 a 8 Gy, 48h) sobre PBMC humanas. O foco principal nao e responder se radiacao "gera Parkinson", mas mapear, com rigor, como a biologia celular e remodelada e em que pontos esse remodelamento se conecta a eixos relevantes para neurodegeneracao.

Objetivos centrais:
1. Caracterizar o efeito dose-resposta global em genes e vias.
2. Identificar assinaturas neurodegenerativas relevantes (sem inferencia diagnostica).
3. Diferenciar vias enriquecidas, empobrecidas e irrelevantes no desenho atual.
4. Fortalecer interpretacao com convergencia entre metodos e robustez entre doadores.

---

## O que este estudo responde (e o que nao responde)

Este estudo responde:
1. Quais sistemas biologicos sobem ou descem com dose de radiacao.
2. Quais achados sao consistentes entre GSVA, FGSEA e CAMERA.
3. Quais sinais se mantem estaveis em analise leave-one-donor-out.
4. Quais padroes sao compatveis com cenario de vulnerabilidade neurodegenerativa.

Este estudo nao responde:
1. Diagnostico de doenca neurodegenerativa.
2. Causalidade clinica individual de Parkinson em humanos.
3. Equivalencia entre assinatura em PBMC e patologia neuronal direta.

---

## Achados principais (mensagem cientifica)

### 1) Radiacao remodela a maquinaria celular de forma ampla e dose-dependente
O efeito nao e localizado em poucos genes isolados. Existe reorganizacao sistemica com ativacao de eixos de dano/inflamacao e queda de eixos de manutencao celular.

### 2) O ponto principal nao e "PD sim ou nao", e sim o padrao de vulnerabilidade biologica
No recorte PD, o sinal e seletivo e nao amplo:
1. Familial PD genes: GSVA logFC por Gy = 0.04548, FDR = 0.01056.
2. PD all: efeito pequeno e nao significativo.

Interpretacao:
1. O dado sugere conexoes neurodegenerativas especificas.
2. O dado nao sustenta assinatura global de doença Parkinson-like em sangue.

### 3) A descoberta mais importante esta no scan global de modulos
Em toda a biblioteca:
1. 236 modulos testados.
2. 75 significativos por GSVA (FDR < 0.05).
3. 63 significativos por FGSEA.
4. 32 significativos por CAMERA.
5. 28 significativos nos tres metodos simultaneamente.

Leitura biologica integrada:
1. Enriquecimento de eixos de estresse, dano e inflamacao.
2. Empobrecimento de eixos de proteostase, ribossomo e manutencao mitocondrial.
3. Perfil altamente compativel com cenario de vulnerabilidade molecular sistemica.

### 4) Escore de risco molecular: tendencia, sem extrapolacao clinica
Score PD-clean:
1. slope = 0.04095 por Gy.
2. p = 0.0743.

Interpretacao:
1. Existe tendencia biologica de risco molecular.
2. Nao ha base para inferencia clinica individual.

### 5) Robustez estatistica e de reproducibilidade sustentam os principais argumentos
1. Concordancia GSVA vs FGSEA: Spearman rho global = 0.908 (n = 236).
2. Concordancia no subconjunto neurodeg-relacionado: rho = 0.904 (n = 134).
3. Modulos com estabilidade LODO alta (consistencia de sinal = 1.0) incluem eixos de ribossomo/proteostase, mtDNA replication e ISR death factors.

---

## Estrutura do dataset

Fonte:
1. NASA GeneLab OSD-157 / GLDS-157.

Desenho:
1. 5 doadores biologicos (a, b, X, Y, Z).
2. 5 doses (0, 0.5, 2, 5, 8 Gy).
3. 25 amostras totais.
4. Endpoint transcriptomico em 48h.

Arquivos principais de entrada:
1. GLDS-157_microarray_*.txt.gz
2. PD2_Mm (1).gmt
3. GLDS-157_array_SampleTable_G4112A_GLmicroarray.csv

---

## Pipeline analitico (fases)

### Fase A: Preprocessamento
1. Leitura Agilent (single-channel).
2. Correcao de background (normexp + offset).
3. Transformacao log2.
4. Normalizacao quantilica.
5. Colapso probe para gene por limma.

### Fase B: Modelagem estatistica
Modelo principal:

Y = beta0 + beta1 * Dose + beta2 * Donor + erro

Componentes:
1. limma para tendencia dose-resposta.
2. duplicateCorrelation para estrutura repetida por doador.
3. arrayWeights para heterogeneidade de arrays.

### Fase C: Enriquecimento e atividade de vias
1. GSVA (atividade por amostra + tendencia por dose).
2. FGSEA (enriquecimento pre-ranqueado).
3. CAMERA (teste competitivo ortogonal).

### Fase D: Score e analise neurodeg-focada
1. Score robusto de risco PD-like sem sobreposicao com nucleo de mitofagia.
2. Teste de tendencia do score com dose.

### Fase E: Scan global de modulos
1. Teste sistematico de todos os modulos analisaveis da biblioteca.
2. Identificacao de vias enriquecidas, empobrecidas e irrelevantes com criterios explicitos.

### Fase F: Robustez e reproducibilidade
1. Efeito com intervalo de confianca.
2. Concordancia GSVA x FGSEA.
3. Distribuicao de p-valores.
4. Leave-one-donor-out (LODO) para estabilidade de sinal.

---

## Evidencias biologicas por classe

### Vias enriquecidas (ativadas)
Criterio: GSVA logFC positivo e FDR abaixo de 0.05.

Exemplos de destaque:
1. ISR / death factors.
2. Antigen presentation.
3. Inflammasome.
4. RAAS tissue damage / PANapoptosis.

Leitura biologica:
1. Aumento de sinalizacao de dano e inflamacoes.
2. Ativacao de resposta adaptativa ao estresse celular.

### Vias empobrecidas (suprimidas)
Criterio: GSVA logFC negativo e FDR abaixo de 0.05.

Exemplos de destaque:
1. Cyto-ribosome assembly/subunits.
2. mtDNA replication.
3. NHEJ (reparo de DNA) em subgrupos.

Leitura biologica:
1. Reducao de manutencao celular e capacidade biossintetica.
2. Maior fragilidade funcional sob estresse persistente.

### Vias irrelevantes no desenho atual
Criterio estrito: efeito muito pequeno + ausencia de suporte convergente entre GSVA/FGSEA/CAMERA.

Interpretacao:
1. Irrelevante aqui significa baixa informatividade neste dataset e nestas condicoes.
2. Nao significa impossibilidade biologica em outros tecidos/tempos.

---

## Resultados por fases (estatisticas e pontuacoes)

### Fase neurodeg-focada
1. Familial PD genes significativo por GSVA.
2. Mito-membrane organization com tendencia positiva, mas sem robustez por multiplicidade.
3. PD all sem efeito robusto.

Arquivos:
1. results/deep_analysis/report_tables/macro_category_impact_summary.csv
2. results/deep_analysis/report_tables/pathway_evidence_table.csv
3. results/deep_analysis/report_tables/risk_score_model_summary.csv

### Fase global
1. Cobertura: 281 parseados, 278 mapeados, 236 testados.
2. Assinatura mista com ativacao imune/estresse e supressao de proteostase/ribossomo.

Arquivos:
1. results/deep_analysis/global_module_scan/module_scan_all.csv
2. results/deep_analysis/global_module_scan/module_scan_by_category.csv
3. results/deep_analysis/global_module_scan/module_scan_top_hits.csv
4. results/deep_analysis/global_module_scan/module_scan_top_abs_effects.csv

### Fase de robustez
1. Concordancia alta entre metodos (rho acima de 0.9).
2. Modulos com estabilidade LODO elevada.

Arquivos:
1. results/deep_analysis/global_module_scan/validated_tables/top_neurodeg_ci95.csv
2. results/deep_analysis/global_module_scan/validated_tables/neurodeg_lodo_stability_summary.csv
3. results/deep_analysis/global_module_scan/validated_tables/neurodeg_leave_one_donor_out_slopes.csv

---

## Conteudo do repositorio

Scripts principais:
1. deep_pd2_analysis_pipeline.R
2. deep_pd2_visual_report.R
3. deep_pd2_ortholog_sensitivity.R
4. deep_pd2_global_module_scan.R
5. deep_pd2_global_extended_figures.R
6. deep_pd2_validated_reproducibility_plots.R

Relatorios:
1. SCIENTIFIC_REPORT_PD2.md
2. EXPLANATORY_EVIDENCE_REPORT_PD2.md
3. Full_Research_Log.md

Pastas de resultados:
1. results/deep_analysis/
2. results/deep_analysis/figures/
3. results/deep_analysis/global_module_scan/
4. results/deep_analysis/global_module_scan/figures/validated/

---

## Como executar

Pipeline principal:
Rscript deep_pd2_analysis_pipeline.R

Figuras neurodeg-focadas:
Rscript deep_pd2_visual_report.R

Sensibilidade ortologa:
Rscript deep_pd2_ortholog_sensitivity.R

Scan global de modulos:
Rscript deep_pd2_global_module_scan.R

Figuras estendidas e robustez:
Rscript deep_pd2_global_extended_figures.R
Rscript deep_pd2_validated_reproducibility_plots.R

---

## Limites do estudo

1. Tamanho amostral de doadores ainda pequeno (n = 5).
2. PBMC e tecido periferico, nao tecido neural.
3. Sobreposicao entre conjuntos de genes pode gerar dependencia entre vias.
4. Mapeamento ortologo online completo depende de conectividade externa.

---

## Conclusao final

A principal contribuicao deste projeto e o mapa sistemico de como a radiacao reorganiza a biologia celular em sangue, com foco em conexoes neurodegenerativas relevantes. O valor cientifico central esta na descricao robusta de vulnerabilidade molecular (eixos de dano/inflamacao versus perda de manutencao celular), e nao em afirmar cenario patologico clinico.

Assim, os dados sustentam um modelo de risco biologico/transcriptomico contextual e reproduzivel, com utilidade para priorizacao mecanistica e estudos de validacao futuros.

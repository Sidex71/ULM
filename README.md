# ULM
An R package for the prediction and identification of multiplets from scRNAseq datasets to infer physical cell-cell interaction networks. Multiplets occur naturally in conventional scRNAseq due to incomplete dissociation during library preparation. They represent cells which are physically connected and interacting in tissues which become sequenced together as they remain unseperated. ULM utilizes a signature-based approach where univariate linear models are fitted over each barcode in a scRNAseq data to assign signature scores. Barcodes are then classified as singlets or multiplets based on their signature scores. Multiplets are those barcodes or cells that are enriched in two or more cell type-specific gene signatures.



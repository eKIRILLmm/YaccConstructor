module YC.GLL.Abstarct.Tests.BioDataPerformance

open VDS.RDF
open VDS.RDF.Parsing

open QuickGraph
open AbstractAnalysis.Common
open Yard.Generators.GLL.AbstractParser
open Yard.Generators.Common.ASTGLL
open Yard.Generators.GLL.ParserCommon

open System.IO
open System.Collections.Generic



let (|Prefix|_|) (p:string) (s:string) =
    if s.StartsWith(p) then Some() else None

let (|Equals|_|) x y = 
    if x = y then Some() else None

let getTokenFromTag tokenizer (tag:string) = 
    match tag with
    | Prefix "Protein_" () -> tokenizer "PROTEIN"
    | Prefix "Gene_" () -> tokenizer "GENE"
    | Prefix "Phenotype_" () -> tokenizer "PHENOTYPE"
    | Equals "belongs_to" () -> tokenizer "BELONGS" 
    | Equals "-belongs_to" () -> tokenizer "RBELONGS" 
    | Equals "codes_for" () -> tokenizer "CODESFOR"
    | Equals "-codes_for" () -> tokenizer "RCODESFOR"
    | Equals "refers_to" () -> tokenizer "REFERS"
    | Prefix "GO_" () -> tokenizer "BIOPROC"
    | Prefix "FamOrDom_" () -> tokenizer "FAM_OR_DOM"
    | Equals "has_FamOrDom" () -> tokenizer "HAS_FAM_OR_DOM"
    | _ -> tokenizer "OTHER"


let loadFromFile (file:string) =
    let g = new Graph()
    if (System.IO.Path.GetExtension file).ToLower() = "ttl"
    then        
        let ttlparser = new TurtleParser()
        ttlparser.Load(g, file)
    else
        FileLoader.Load(g, file)       
    g

let vMap = new System.Collections.Generic.Dictionary<_,_>()
let mutable idV = -1
let getId v = if vMap.ContainsKey v then vMap.[v] else (idV <- idV + 2; vMap.Add(v, idV); idV)
 

let getEdgesFromHomologene file =
    let lines = File.ReadAllLines(file)
    [|
        for l in lines->
            let elems = l.Split('\t')
            let homoloGeneGroup = "HomoloGene_" + elems.[0] 
            let gene = "Gene_" + elems.[2]
            [|
                getId gene + 1, "is_homologous_to", getId homoloGeneGroup;
                getId homoloGeneGroup, homoloGeneGroup, getId homoloGeneGroup + 1;
            |]
    |] |> Array.concat

let getEdgesFromKeggPath file = 
    let lines = File.ReadAllLines(file)
    [|
        for l in lines->
            let elems = l.Split('\t')
            let pathw = "Pathway_" + elems.[0]
            let namePath = "name_pathway" + elems.[1]
            [|
                getId pathw, pathw, getId pathw + 1
                getId pathw + 1, "name path", getId namePath;
                getId namePath, namePath, getId namePath + 1;
            |]
    |] |> Array.concat

let getEdgesFromKeggGeneMap file =
    let lines = File.ReadAllLines(file)
    [|
        for l in lines->
            let elems = l.Split('\t')
            let gene = "Gene_" + elems.[0]
            let pathw = "Pathway_" + elems.[1]
            [|
                getId gene + 1, "has_Pathway", getId pathw;
            |]
    |] |> Array.concat

let getEdgesFromInterpro file =
    let lines = File.ReadAllLines(file)
    [|
        for l in lines->
            let elems = l.Split('\t')
            let ipro = "FamOrDom_" + elems.[0];
            let protein_count = elems.[1]
            let short_name = elems.[2]
            let iproType = elems.[3]
            let name = elems.[4]
            let iproId = getId ipro
            [|
                iproId, ipro, iproId + 1;

                iproId + 1, "potein_count", getId protein_count;
                getId protein_count, protein_count, getId protein_count + 1;

                iproId + 1, "iprot_short_name", getId short_name;
                getId short_name, short_name, getId short_name + 1;

                iproId + 1, "iprot_type", getId iproType;
                getId iproType, iproType, getId iproType + 1;

                iproId + 1, "iprot_name", getId name;
                getId name, name, getId name + 1;
            |]
    |] |>Array.concat
    

let getEdgesFromEntrezGene file =
    let lines = File.ReadAllLines(file)
    [|
        for l in lines ->
            let elems = l.Split('\t')

            let tax_id = elems.[0]
            let Gene = "Gene_" + elems.[1]
            let Symbol = elems.[2]
            let LocusTag = elems.[3]
            let Synonyms = elems.[4]
            let dbXrefs = elems.[5]
            let chromosome = elems.[6]
            let map_location = elems.[7]
            let description = elems.[8]
            let type_of_gene = elems.[9]
            let Symbol_from_nomenclature_authority = elems.[10]
            let Full_name_from_nomenclature_authority = elems.[11]
            let Nomenclature_status = elems.[12]
            let Other_designations = elems.[13]
            let Modification_date = elems.[14]

            let geneId = getId Gene
            //printfn "%A, %A, %A" geneId Gene (geneId+1)
            [|
                geneId, Gene, geneId + 1;

                geneId + 1, "tax_id", getId tax_id;
                getId tax_id, tax_id,  (getId tax_id) + 1;

                geneId + 1, "Symbol", getId Symbol;
                getId Symbol, Symbol,  (getId Symbol) + 1;

                geneId + 1, "LocusTag", getId LocusTag;
                getId LocusTag, LocusTag,  (getId LocusTag) + 1;

                geneId + 1, "Synonyms", getId Synonyms;
                getId Synonyms, Synonyms,  (getId Synonyms) + 1;

                geneId + 1, "dbXrefs", getId dbXrefs;
                getId dbXrefs, dbXrefs,  (getId dbXrefs) + 1;

                geneId + 1, "chromosome", getId chromosome;
                getId chromosome, chromosome,  (getId chromosome) + 1;

                geneId + 1, "description", getId description;
                getId description, description,  (getId description) + 1;

                geneId + 1, "type_of_gene", getId type_of_gene;
                getId type_of_gene, type_of_gene,  (getId type_of_gene) + 1;

                geneId + 1, "map_location", getId map_location;
                getId map_location, map_location,  (getId map_location) + 1;

                geneId + 1, "Symbol_from_nomenclature_authority", getId Symbol_from_nomenclature_authority;
                getId Symbol_from_nomenclature_authority, Symbol_from_nomenclature_authority,  (getId Symbol_from_nomenclature_authority) + 1;

                geneId + 1, "Full_name_from_nomenclature_authority", getId Full_name_from_nomenclature_authority;
                getId Full_name_from_nomenclature_authority, Full_name_from_nomenclature_authority,  (getId Full_name_from_nomenclature_authority) + 1;

                geneId + 1, "Nomenclature_status", getId Nomenclature_status;
                getId Nomenclature_status, Nomenclature_status,  (getId Nomenclature_status) + 1;

                geneId + 1, "Other_designations", getId Other_designations;
                getId Other_designations, Other_designations,  (getId Other_designations) + 1;

                geneId + 1, "Modification_date", getId Modification_date;
                getId Modification_date, Modification_date,  (getId Modification_date) + 1;

            |]
    |] |> Array.concat
    

let getEdgesFromGO file =
    let g = loadFromFile file

    let edg (f: VDS.RDF.INode) (t: VDS.RDF.INode) (l: VDS.RDF.INode) = 
        match f, t, l with
        | f, t, l when f.ToString().StartsWith "http://purl.obolibrary.org/obo" ->
            let GOid = (f :?> UriNode).Uri.Segments.[(f :?> UriNode).Uri.Segments.Length - 1]
            [|
                getId GOid, GOid, getId GOid + 1;
                getId GOid + 1, l.ToString(), getId (t.ToString());
                getId (t.ToString()), t.ToString(), getId (t.ToString()) + 1;
            |]
        | _ -> 
            let fId = getId (f.ToString())
            let tId = getId (t.ToString())
            [|fId, f.ToString(), fId + 1;
            fId + 1, l.ToString(), tId;
            tId, t.ToString(), tId + 1|]

    let edgs = [|for t in g.Triples -> edg t.Object t.Subject t.Predicate |] |> Array.concat
    edgs

let getEdgesFromUniprot file =
    let g = loadFromFile file
    let triples = g.Triples.Count

    let mutable curProt = ""
    let mutable curGene = ""
    let edg (f: VDS.RDF.INode) (t: VDS.RDF.INode) (l: VDS.RDF.INode) = 
        match f, t, l with
        | f, t, l when f.ToString() = "http://purl.uniprot.org/core/Protein" -> 
            curProt <- "Protein_" + (t :?> UriNode).Uri.Segments.[(t :?> UriNode).Uri.Segments.Length - 1]
            let fId = getId (f.ToString())
            let tId = getId (t.ToString())
            let protId = getId curProt
            //printfn "%A, %A, %A" protId curProt (protId + 1)
            [|
                fId, f.ToString(), fId + 1;
                fId + 1, l.ToString(), protId;
                protId, curProt, protId + 1
            |]
        | f, t, l when f.ToString().StartsWith "http://purl.uniprot.org/interpro" -> 
            let interpro = "FamOrDom_" + (f :?> UriNode).Uri.Segments.[(f :?> UriNode).Uri.Segments.Length - 1]
            let interproID = getId interpro
            let protId = getId curProt
            //printfn "%A, %A, %A" fId curProt tId
            [|
                protId + 1, "has_FamOrDom", interproID;
                interproID, interpro, interproID + 1
            |]
            
        | f, t, l when l.ToString() = "http://purl.uniprot.org/core/classifiedWith" && (f :?> UriNode).Uri.Segments.[(f :?> UriNode).Uri.Segments.Length - 1].StartsWith("GO") ->
            let GOid = (f :?> UriNode).Uri.Segments.[(f :?> UriNode).Uri.Segments.Length - 1]
            let fId = getId curProt
            let tId = getId GOid
            //printfn "%A, %A, %A" fId GOid tId

            [|
                fId + 1, "belongs_to", tId;
                tId + 1, "-belongs_to", fId;
            |]

        | f, t, l when f.ToString() = "http://purl.uniprot.org/database/GeneID" ->
            let gene = "Gene_" + (t :?> UriNode).Uri.Segments.[(t :?> UriNode).Uri.Segments.Length - 1]
            curGene <- gene
            let protId = getId curProt
            let geneId = getId gene
            //printfn "%A, %A, %A" geneId gene (geneId+1)
            [|
                protId + 1, "-codes_for", geneId;
                geneId + 1, "codes_for", protId;
            |]

        | f, t, l when f.ToString() = "phenotype" -> 
            let MIMid = "Phenotype_" + (t :?> UriNode).Uri.Segments.[(t :?> UriNode).Uri.Segments.Length - 1]
            let fId = getId curGene
            let tId = getId MIMid
            //printfn "%A, %A, %A" fId MIMid tId
            [|
                fId + 1, "refers_to", tId;
                tId, MIMid, tId + 1;
            |]

        | _ -> let fId = getId (f.ToString())
               let tId = getId (t.ToString())
               [|
                   fId, f.ToString(), fId + 1;
                   fId + 1, l.ToString(), tId;
                   tId, t.ToString(), tId + 1
               |]
               
    
    let edgs = [|for t in g.Triples -> edg t.Object t.Subject t.Predicate |] |> Array.concat
    edgs

let getParseInputGraph files =
//    let getEdges f = 
//        match System.IO.Path.GetFileName f with
//        | "go.owl" -> getEdgesFromGO f
//        | "uniprot-for_test.rdf" 
//        | "uniprot-homo_sapiens.rdf" 
//        | "uniprot-c_elegans.rdf" 
//        | "uniprot-melanogaster.rdf" 
//        | "uniprot-mus_musculus.rdf" 
//        | "uniprot-rattus_norvegicus.rdf" -> getEdgesFromUniprot f
//        | "Caenorhabditis_elegans.gene_info"
//        | "Drosophila_melanogaster.gene_info"
//        | "Mus_musculus.gene_info"
//        | "Rattus_norvegicus.gene_info"
//        | "Homo_sapiens.gene_info" -> getEdgesFromEntrezGene f
//        | "interpro.txt" -> getEdgesFromInterpro f
//        | "pathways.keg" -> getEdgesFromKeggPath f
//        | "geneToPath.txt" -> getEdgesFromKeggGeneMap f
//        | "homologene.data.txt" -> getEdgesFromHomologene f
//        | _ -> [||]
//
//    let edgs = files 
//            |> Array.map getEdges 
//            |> Array.concat 
    
    let edgs = 
        [|
            1, "Gene_1", 2;
            2, "codes_for", 3;
            3, "Protein_1", 4;
            4, "belongs_to", 5;
            5, "GO_1", 6;
            6, "-belongs_to", 7;
            7, "Protein_2", 8;
            8, "-codes_for", 9;
            9, "Gene_2", 10;

//            4, "-codes_for", 1;
//            6, "-belongs_to", 3;
//            8, "belongs_to", 5;
//            10, "codes_for", 7;
        |]

    let allVs = edgs |> Array.collect (fun (f,l,t) -> [|f * 1<positionInInput>; t * 1<positionInInput>|]) |> Set.ofArray |> Array.ofSeq
    let eofV = allVs.Length
        
    let graph = new SimpleInputGraph<_>(allVs, getTokenFromTag (fun x -> (int) GLL.GPPerf1.stringToToken.[x]))
    
    edgs
    |> Array.collect (fun (f,l,t) -> [|new ParserEdge<_>(f, t, l)|])
    |> graph.AddVerticesAndEdgeRange
    |> ignore



    graph, graph.EdgeCount
        
let processFiles files =
    let g1, edges = 
        getParseInputGraph files 
//    for e in g1.Edges do
//        printfn "%A; %A" (e.Tag) (e.ToString())
    let cnt = 1
    let start = System.DateTime.Now
    printfn "%A" edges
    let root1 =
        Yard.Generators.GLL.AbstractParser.getAllRangesForStartState GLL.GPPerf1.parserSource g1
        |> Seq.length
    
    let time1 = (System.DateTime.Now - start).TotalMilliseconds / (float cnt)


    edges, time1, root1

let performTests () =
    let basePath = @"..\..\..\data\BioData"
    let files = System.IO.Directory.GetFiles basePath    
    files 
    |> processFiles
    |> printfn "%A"
// Implementation file for parser generated by fsyacc
module AntlrToYard.Parser
#nowarn "64";; // turn off warnings that type variables used in production annotations are instantiated to concrete type
open Yard.Core.IL
open Microsoft.FSharp.Text.Lexing
open Microsoft.FSharp.Text.Parsing.ParseHelpers
# 1 "Parser.fsy"


open Yard.Core.IL.Production

(* Run with fsyacc.exe --module AntlrToYard.Parser --open Yard.Core.IL Parser.fsy *)


(*
Expr: ID { Val($1) }
     | INT {  Int($1)  }
     | FLOAT {  Float($1)  }
     | DECR LPAREN Expr RPAREN {  Decr($3)  }


 Stmt: ID ASSIGN Expr { Assign($1,$3) }
     | WHILE Expr DO Stmt { While($2,$4) }
     | BEGIN StmtList END { Seq(List.rev($2)) }
     | IF Expr THEN Stmt { IfThen($2,$4) }
     | IF Expr THEN Stmt ELSE Stmt { IfThenElse($2,$4,$6) }
     | PRINT Expr { Print($2) }


 StmtList: Stmt { [$1] }
        | StmtList SEMI Stmt { $3 :: $1  }

*)

let makeModifiedRule innerProduction modifier =
    match modifier with
    | "+" -> PSome(innerProduction)
    | "*" -> PMany(innerProduction)
    | "?" -> POpt(innerProduction)
    | "!" -> innerProduction // Not included in AST
    | "" -> innerProduction

let makePSeq (productionList, actionCode) =
    PSeq( List.map (fun prod -> {omit = false; rule = prod; binding = None; checker = None;}) productionList , actionCode )

let termCount = ref 0
let generateNewName =
    termCount := !termCount + 1
    sprintf "TERMINAL_%i" !termCount 
let terminals = new System.Collections.Generic.Dictionary<string, string>()
let makeToken (identifier, pos) descr =
    let newName = if identifier="" then generateNewName else identifier
    terminals.[newName] <- descr
    PToken(identifier, pos)
        

# 57 "Parser.fs"
// This type is the type of tokens accepted by the parser
type token = 
  | DOUBLE_DOT
  | TILDE
  | EXCLAMATION
  | QUESTION
  | SEMICOLON
  | COLON
  | PLUS
  | STAR
  | EQUAL
  | BAR
  | RPAREN
  | LPAREN
  | TERMINAL of (Source.t)
  | LITERAL of (Source.t)
  | IDENTIFIER of (Source.t)
  | T_OPTIONS
  | T_GRAMMAR
  | EOF
  | ACTION_CODE of (Source.t)
  | ACTION_NAME of (Source.t)
  | SCOPE_NAME of (Source.t)
  | SINGLELINE_COMMENT of (Source.t)
  | MULTILINE_COMMENT of (Source.t)
// This type is used to give symbolic names to token indexes, useful for error messages
type tokenId = 
    | TOKEN_DOUBLE_DOT
    | TOKEN_TILDE
    | TOKEN_EXCLAMATION
    | TOKEN_QUESTION
    | TOKEN_SEMICOLON
    | TOKEN_COLON
    | TOKEN_PLUS
    | TOKEN_STAR
    | TOKEN_EQUAL
    | TOKEN_BAR
    | TOKEN_RPAREN
    | TOKEN_LPAREN
    | TOKEN_TERMINAL
    | TOKEN_LITERAL
    | TOKEN_IDENTIFIER
    | TOKEN_T_OPTIONS
    | TOKEN_T_GRAMMAR
    | TOKEN_EOF
    | TOKEN_ACTION_CODE
    | TOKEN_ACTION_NAME
    | TOKEN_SCOPE_NAME
    | TOKEN_SINGLELINE_COMMENT
    | TOKEN_MULTILINE_COMMENT
    | TOKEN_end_of_input
    | TOKEN_error
// This type is used to give symbolic names to token indexes, useful for error messages
type nonTerminalId = 
    | NONTERM__startParseAntlr
    | NONTERM_ParseAntlr
    | NONTERM_TopLevelDefs
    | NONTERM_TopLevelDef
    | NONTERM_Rule
    | NONTERM_TerminalRule
    | NONTERM_Options
    | NONTERM_RuleBody
    | NONTERM_Alt
    | NONTERM_ActionCodeOptional
    | NONTERM_Seq
    | NONTERM_Modifier
    | NONTERM_SimpleProduction
    | NONTERM_RuleString
    | NONTERM_RulePart

// This function maps tokens to integers indexes
let tagOfToken (t:token) = 
  match t with
  | DOUBLE_DOT  -> 0 
  | TILDE  -> 1 
  | EXCLAMATION  -> 2 
  | QUESTION  -> 3 
  | SEMICOLON  -> 4 
  | COLON  -> 5 
  | PLUS  -> 6 
  | STAR  -> 7 
  | EQUAL  -> 8 
  | BAR  -> 9 
  | RPAREN  -> 10 
  | LPAREN  -> 11 
  | TERMINAL _ -> 12 
  | LITERAL _ -> 13 
  | IDENTIFIER _ -> 14 
  | T_OPTIONS  -> 15 
  | T_GRAMMAR  -> 16 
  | EOF  -> 17 
  | ACTION_CODE _ -> 18 
  | ACTION_NAME _ -> 19 
  | SCOPE_NAME _ -> 20 
  | SINGLELINE_COMMENT _ -> 21 
  | MULTILINE_COMMENT _ -> 22 

// This function maps integers indexes to symbolic token ids
let tokenTagToTokenId (tokenIdx:int) = 
  match tokenIdx with
  | 0 -> TOKEN_DOUBLE_DOT 
  | 1 -> TOKEN_TILDE 
  | 2 -> TOKEN_EXCLAMATION 
  | 3 -> TOKEN_QUESTION 
  | 4 -> TOKEN_SEMICOLON 
  | 5 -> TOKEN_COLON 
  | 6 -> TOKEN_PLUS 
  | 7 -> TOKEN_STAR 
  | 8 -> TOKEN_EQUAL 
  | 9 -> TOKEN_BAR 
  | 10 -> TOKEN_RPAREN 
  | 11 -> TOKEN_LPAREN 
  | 12 -> TOKEN_TERMINAL 
  | 13 -> TOKEN_LITERAL 
  | 14 -> TOKEN_IDENTIFIER 
  | 15 -> TOKEN_T_OPTIONS 
  | 16 -> TOKEN_T_GRAMMAR 
  | 17 -> TOKEN_EOF 
  | 18 -> TOKEN_ACTION_CODE 
  | 19 -> TOKEN_ACTION_NAME 
  | 20 -> TOKEN_SCOPE_NAME 
  | 21 -> TOKEN_SINGLELINE_COMMENT 
  | 22 -> TOKEN_MULTILINE_COMMENT 
  | 25 -> TOKEN_end_of_input
  | 23 -> TOKEN_error
  | _ -> failwith "tokenTagToTokenId: bad token"

/// This function maps production indexes returned in syntax errors to strings representing the non terminal that would be produced by that production
let prodIdxToNonTerminal (prodIdx:int) = 
  match prodIdx with
    | 0 -> NONTERM__startParseAntlr 
    | 1 -> NONTERM_ParseAntlr 
    | 2 -> NONTERM_TopLevelDefs 
    | 3 -> NONTERM_TopLevelDefs 
    | 4 -> NONTERM_TopLevelDef 
    | 5 -> NONTERM_TopLevelDef 
    | 6 -> NONTERM_TopLevelDef 
    | 7 -> NONTERM_Rule 
    | 8 -> NONTERM_TerminalRule 
    | 9 -> NONTERM_Options 
    | 10 -> NONTERM_Options 
    | 11 -> NONTERM_RuleBody 
    | 12 -> NONTERM_RuleBody 
    | 13 -> NONTERM_Alt 
    | 14 -> NONTERM_Alt 
    | 15 -> NONTERM_ActionCodeOptional 
    | 16 -> NONTERM_ActionCodeOptional 
    | 17 -> NONTERM_Seq 
    | 18 -> NONTERM_Seq 
    | 19 -> NONTERM_Modifier 
    | 20 -> NONTERM_Modifier 
    | 21 -> NONTERM_Modifier 
    | 22 -> NONTERM_Modifier 
    | 23 -> NONTERM_Modifier 
    | 24 -> NONTERM_SimpleProduction 
    | 25 -> NONTERM_SimpleProduction 
    | 26 -> NONTERM_SimpleProduction 
    | 27 -> NONTERM_SimpleProduction 
    | 28 -> NONTERM_RuleString 
    | 29 -> NONTERM_RuleString 
    | 30 -> NONTERM_RulePart 
    | 31 -> NONTERM_RulePart 
    | 32 -> NONTERM_RulePart 
    | 33 -> NONTERM_RulePart 
    | 34 -> NONTERM_RulePart 
    | 35 -> NONTERM_RulePart 
    | 36 -> NONTERM_RulePart 
    | 37 -> NONTERM_RulePart 
    | 38 -> NONTERM_RulePart 
    | 39 -> NONTERM_RulePart 
    | 40 -> NONTERM_RulePart 
    | 41 -> NONTERM_RulePart 
    | 42 -> NONTERM_RulePart 
    | 43 -> NONTERM_RulePart 
    | 44 -> NONTERM_RulePart 
    | 45 -> NONTERM_RulePart 
    | _ -> failwith "prodIdxToNonTerminal: bad production index"

let _fsyacc_endOfInputTag = 25 
let _fsyacc_tagOfErrorTerminal = 23

// This function gets the name of a token as a string
let token_to_string (t:token) = 
  match t with 
  | DOUBLE_DOT  -> "DOUBLE_DOT" 
  | TILDE  -> "TILDE" 
  | EXCLAMATION  -> "EXCLAMATION" 
  | QUESTION  -> "QUESTION" 
  | SEMICOLON  -> "SEMICOLON" 
  | COLON  -> "COLON" 
  | PLUS  -> "PLUS" 
  | STAR  -> "STAR" 
  | EQUAL  -> "EQUAL" 
  | BAR  -> "BAR" 
  | RPAREN  -> "RPAREN" 
  | LPAREN  -> "LPAREN" 
  | TERMINAL _ -> "TERMINAL" 
  | LITERAL _ -> "LITERAL" 
  | IDENTIFIER _ -> "IDENTIFIER" 
  | T_OPTIONS  -> "T_OPTIONS" 
  | T_GRAMMAR  -> "T_GRAMMAR" 
  | EOF  -> "EOF" 
  | ACTION_CODE _ -> "ACTION_CODE" 
  | ACTION_NAME _ -> "ACTION_NAME" 
  | SCOPE_NAME _ -> "SCOPE_NAME" 
  | SINGLELINE_COMMENT _ -> "SINGLELINE_COMMENT" 
  | MULTILINE_COMMENT _ -> "MULTILINE_COMMENT" 

// This function gets the data carried by a token as an object
let _fsyacc_dataOfToken (t:token) = 
  match t with 
  | DOUBLE_DOT  -> (null : System.Object) 
  | TILDE  -> (null : System.Object) 
  | EXCLAMATION  -> (null : System.Object) 
  | QUESTION  -> (null : System.Object) 
  | SEMICOLON  -> (null : System.Object) 
  | COLON  -> (null : System.Object) 
  | PLUS  -> (null : System.Object) 
  | STAR  -> (null : System.Object) 
  | EQUAL  -> (null : System.Object) 
  | BAR  -> (null : System.Object) 
  | RPAREN  -> (null : System.Object) 
  | LPAREN  -> (null : System.Object) 
  | TERMINAL _fsyacc_x -> Microsoft.FSharp.Core.Operators.box _fsyacc_x 
  | LITERAL _fsyacc_x -> Microsoft.FSharp.Core.Operators.box _fsyacc_x 
  | IDENTIFIER _fsyacc_x -> Microsoft.FSharp.Core.Operators.box _fsyacc_x 
  | T_OPTIONS  -> (null : System.Object) 
  | T_GRAMMAR  -> (null : System.Object) 
  | EOF  -> (null : System.Object) 
  | ACTION_CODE _fsyacc_x -> Microsoft.FSharp.Core.Operators.box _fsyacc_x 
  | ACTION_NAME _fsyacc_x -> Microsoft.FSharp.Core.Operators.box _fsyacc_x 
  | SCOPE_NAME _fsyacc_x -> Microsoft.FSharp.Core.Operators.box _fsyacc_x 
  | SINGLELINE_COMMENT _fsyacc_x -> Microsoft.FSharp.Core.Operators.box _fsyacc_x 
  | MULTILINE_COMMENT _fsyacc_x -> Microsoft.FSharp.Core.Operators.box _fsyacc_x 
let _fsyacc_gotos = [| 0us; 65535us; 1us; 65535us; 0us; 1us; 1us; 65535us; 0us; 2us; 2us; 65535us; 0us; 4us; 2us; 5us; 2us; 65535us; 0us; 6us; 2us; 6us; 2us; 65535us; 0us; 7us; 2us; 7us; 2us; 65535us; 10us; 11us; 15us; 16us; 2us; 65535us; 12us; 13us; 30us; 23us; 4us; 65535us; 12us; 22us; 24us; 25us; 26us; 28us; 30us; 22us; 1us; 65535us; 26us; 27us; 4us; 65535us; 12us; 26us; 24us; 26us; 26us; 26us; 30us; 26us; 2us; 65535us; 31us; 32us; 33us; 34us; 4us; 65535us; 12us; 33us; 24us; 33us; 26us; 33us; 30us; 33us; 2us; 65535us; 17us; 18us; 44us; 45us; 2us; 65535us; 17us; 44us; 44us; 44us; |]
let _fsyacc_sparseGotoTableRowOffsets = [|0us; 1us; 3us; 5us; 8us; 11us; 14us; 17us; 20us; 25us; 27us; 32us; 35us; 40us; 43us; |]
let _fsyacc_stateToProdIdxsTableElements = [| 1us; 0us; 1us; 0us; 2us; 1us; 3us; 1us; 1us; 1us; 2us; 1us; 3us; 1us; 4us; 1us; 5us; 1us; 6us; 1us; 6us; 1us; 7us; 1us; 7us; 1us; 7us; 2us; 7us; 12us; 1us; 7us; 1us; 8us; 1us; 8us; 1us; 8us; 1us; 8us; 1us; 8us; 1us; 10us; 1us; 10us; 1us; 11us; 2us; 12us; 17us; 1us; 12us; 1us; 12us; 2us; 13us; 14us; 1us; 13us; 1us; 14us; 1us; 16us; 1us; 17us; 1us; 17us; 1us; 17us; 1us; 18us; 1us; 18us; 1us; 19us; 1us; 20us; 1us; 21us; 1us; 22us; 1us; 24us; 2us; 25us; 26us; 1us; 25us; 1us; 25us; 1us; 27us; 2us; 28us; 29us; 1us; 29us; 1us; 30us; 1us; 31us; 1us; 32us; 1us; 33us; 1us; 34us; 1us; 35us; 1us; 36us; 1us; 37us; 1us; 38us; 1us; 39us; 1us; 40us; 1us; 41us; 1us; 42us; 1us; 43us; 1us; 44us; 1us; 45us; |]
let _fsyacc_stateToProdIdxsTableRowOffsets = [|0us; 2us; 4us; 7us; 9us; 11us; 13us; 15us; 17us; 19us; 21us; 23us; 25us; 27us; 30us; 32us; 34us; 36us; 38us; 40us; 42us; 44us; 46us; 48us; 51us; 53us; 55us; 58us; 60us; 62us; 64us; 66us; 68us; 70us; 72us; 74us; 76us; 78us; 80us; 82us; 84us; 87us; 89us; 91us; 93us; 96us; 98us; 100us; 102us; 104us; 106us; 108us; 110us; 112us; 114us; 116us; 118us; 120us; 122us; 124us; 126us; 128us; |]
let _fsyacc_action_rows = 62
let _fsyacc_actionTableElements = [|3us; 32768us; 12us; 15us; 14us; 10us; 15us; 8us; 0us; 49152us; 4us; 32768us; 12us; 15us; 14us; 10us; 15us; 8us; 17us; 3us; 0us; 16385us; 0us; 16386us; 0us; 16387us; 0us; 16388us; 0us; 16389us; 1us; 32768us; 18us; 9us; 0us; 16390us; 1us; 16393us; 15us; 20us; 1us; 32768us; 5us; 12us; 4us; 32768us; 11us; 30us; 12us; 43us; 13us; 40us; 14us; 39us; 2us; 32768us; 4us; 14us; 9us; 24us; 0us; 16391us; 1us; 16393us; 15us; 20us; 1us; 32768us; 5us; 17us; 16us; 32768us; 0us; 46us; 1us; 47us; 2us; 49us; 3us; 48us; 6us; 50us; 7us; 51us; 8us; 52us; 9us; 53us; 10us; 54us; 11us; 55us; 12us; 58us; 13us; 56us; 14us; 57us; 18us; 59us; 19us; 60us; 20us; 61us; 1us; 32768us; 4us; 19us; 0us; 16392us; 1us; 32768us; 18us; 21us; 0us; 16394us; 0us; 16395us; 2us; 32768us; 9us; 24us; 10us; 31us; 4us; 32768us; 11us; 30us; 12us; 43us; 13us; 40us; 14us; 39us; 0us; 16396us; 5us; 16399us; 11us; 30us; 12us; 43us; 13us; 40us; 14us; 39us; 18us; 29us; 0us; 16397us; 0us; 16398us; 0us; 16400us; 4us; 32768us; 11us; 30us; 12us; 43us; 13us; 40us; 14us; 39us; 4us; 16407us; 2us; 38us; 3us; 37us; 6us; 35us; 7us; 36us; 0us; 16401us; 4us; 16407us; 2us; 38us; 3us; 37us; 6us; 35us; 7us; 36us; 0us; 16402us; 0us; 16403us; 0us; 16404us; 0us; 16405us; 0us; 16406us; 0us; 16408us; 1us; 16410us; 0us; 41us; 1us; 32768us; 13us; 42us; 0us; 16409us; 0us; 16411us; 16us; 16412us; 0us; 46us; 1us; 47us; 2us; 49us; 3us; 48us; 6us; 50us; 7us; 51us; 8us; 52us; 9us; 53us; 10us; 54us; 11us; 55us; 12us; 58us; 13us; 56us; 14us; 57us; 18us; 59us; 19us; 60us; 20us; 61us; 0us; 16413us; 0us; 16414us; 0us; 16415us; 0us; 16416us; 0us; 16417us; 0us; 16418us; 0us; 16419us; 0us; 16420us; 0us; 16421us; 0us; 16422us; 0us; 16423us; 0us; 16424us; 0us; 16425us; 0us; 16426us; 0us; 16427us; 0us; 16428us; 0us; 16429us; |]
let _fsyacc_actionTableRowOffsets = [|0us; 4us; 5us; 10us; 11us; 12us; 13us; 14us; 15us; 17us; 18us; 20us; 22us; 27us; 30us; 31us; 33us; 35us; 52us; 54us; 55us; 57us; 58us; 59us; 62us; 67us; 68us; 74us; 75us; 76us; 77us; 82us; 87us; 88us; 93us; 94us; 95us; 96us; 97us; 98us; 99us; 101us; 103us; 104us; 105us; 122us; 123us; 124us; 125us; 126us; 127us; 128us; 129us; 130us; 131us; 132us; 133us; 134us; 135us; 136us; 137us; 138us; |]
let _fsyacc_reductionSymbolCounts = [|1us; 2us; 1us; 2us; 1us; 1us; 2us; 5us; 5us; 0us; 2us; 1us; 3us; 2us; 2us; 0us; 1us; 4us; 2us; 1us; 1us; 1us; 1us; 0us; 1us; 3us; 1us; 1us; 1us; 2us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; 1us; |]
let _fsyacc_productionToNonTerminalTable = [|0us; 1us; 2us; 2us; 3us; 3us; 3us; 4us; 5us; 6us; 6us; 7us; 7us; 8us; 8us; 9us; 9us; 10us; 10us; 11us; 11us; 11us; 11us; 11us; 12us; 12us; 12us; 12us; 13us; 13us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; 14us; |]
let _fsyacc_immediateActions = [|65535us; 49152us; 65535us; 16385us; 16386us; 16387us; 16388us; 16389us; 65535us; 16390us; 65535us; 65535us; 65535us; 65535us; 16391us; 65535us; 65535us; 65535us; 65535us; 16392us; 65535us; 16394us; 16395us; 65535us; 65535us; 16396us; 65535us; 16397us; 16398us; 16400us; 65535us; 65535us; 16401us; 65535us; 16402us; 16403us; 16404us; 16405us; 16406us; 16408us; 65535us; 65535us; 16409us; 16411us; 65535us; 16413us; 16414us; 16415us; 16416us; 16417us; 16418us; 16419us; 16420us; 16421us; 16422us; 16423us; 16424us; 16425us; 16426us; 16427us; 16428us; 16429us; |]
let _fsyacc_reductions ()  =    [| 
# 303 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : (Source.t, Source.t)Grammar.t * (string, string)System.Collections.Generic.Dictionary)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
                      raise (Microsoft.FSharp.Text.Parsing.Accept(Microsoft.FSharp.Core.Operators.box _1))
                   )
                 : '_startParseAntlr));
# 312 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'TopLevelDefs)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 83 "Parser.fsy"
                                                    ((_1), terminals)  
                   )
# 83 "Parser.fsy"
                 : (Source.t, Source.t)Grammar.t * (string, string)System.Collections.Generic.Dictionary));
# 323 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'TopLevelDef)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 85 "Parser.fsy"
                                                 _1 
                   )
# 85 "Parser.fsy"
                 : 'TopLevelDefs));
# 334 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'TopLevelDefs)) in
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : 'TopLevelDef)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 86 "Parser.fsy"
                                                      _1 @ _2 
                   )
# 86 "Parser.fsy"
                 : 'TopLevelDefs));
# 346 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'Rule)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 88 "Parser.fsy"
                                         [_1] 
                   )
# 88 "Parser.fsy"
                 : 'TopLevelDef));
# 357 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'TerminalRule)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 89 "Parser.fsy"
                                          [] 
                   )
# 89 "Parser.fsy"
                 : 'TopLevelDef));
# 368 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 90 "Parser.fsy"
                                                   [] 
                   )
# 90 "Parser.fsy"
                 : 'TopLevelDef));
# 379 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : 'Options)) in
            let _4 = (let data = parseState.GetInput(4) in (Microsoft.FSharp.Core.Operators.unbox data : 'RuleBody)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 92 "Parser.fsy"
                                                                         { new Rule.t<Source.t, Source.t> with name = fst(_1) and args = [] and body = _4 and _public = false and metaArgs = [] } 
                   )
# 92 "Parser.fsy"
                 : 'Rule));
# 392 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : 'Options)) in
            let _4 = (let data = parseState.GetInput(4) in (Microsoft.FSharp.Core.Operators.unbox data : 'RuleString)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 94 "Parser.fsy"
                                                                                 makeToken (_1) (List.fold (fun acc elem -> acc+" "+elem) "" (_4)); 
                   )
# 94 "Parser.fsy"
                 : 'TerminalRule));
# 405 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 96 "Parser.fsy"
                               
                   )
# 96 "Parser.fsy"
                 : 'Options));
# 415 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 96 "Parser.fsy"
                                                          
                   )
# 96 "Parser.fsy"
                 : 'Options));
# 426 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'Alt)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 98 "Parser.fsy"
                                     makePSeq (_1) 
                   )
# 98 "Parser.fsy"
                 : 'RuleBody));
# 437 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'RuleBody)) in
            let _3 = (let data = parseState.GetInput(3) in (Microsoft.FSharp.Core.Operators.unbox data : 'Alt)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 98 "Parser.fsy"
                                                                          PAlt(_1, makePSeq (_3)) 
                   )
# 98 "Parser.fsy"
                 : 'RuleBody));
# 449 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'Seq)) in
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : 'ActionCodeOptional)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 100 "Parser.fsy"
                                                   ([_1], _2) 
                   )
# 100 "Parser.fsy"
                 : 'Alt));
# 461 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'Seq)) in
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : 'Alt)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 100 "Parser.fsy"
                                                                            _1 :: fst(_2), None 
                   )
# 100 "Parser.fsy"
                 : 'Alt));
# 473 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 102 "Parser.fsy"
                                           None 
                   )
# 102 "Parser.fsy"
                 : 'ActionCodeOptional));
# 483 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 103 "Parser.fsy"
                                         Some(_1) 
                   )
# 103 "Parser.fsy"
                 : 'ActionCodeOptional));
# 494 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : 'RuleBody)) in
            let _4 = (let data = parseState.GetInput(4) in (Microsoft.FSharp.Core.Operators.unbox data : 'Modifier)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 105 "Parser.fsy"
                                                            makeModifiedRule _2 _4 
                   )
# 105 "Parser.fsy"
                 : 'Seq));
# 506 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'SimpleProduction)) in
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : 'Modifier)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 106 "Parser.fsy"
                                                       makeModifiedRule (_1) _2 
                   )
# 106 "Parser.fsy"
                 : 'Seq));
# 518 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 108 "Parser.fsy"
                                      "+" 
                   )
# 108 "Parser.fsy"
                 : 'Modifier));
# 528 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 109 "Parser.fsy"
                                  "*" 
                   )
# 109 "Parser.fsy"
                 : 'Modifier));
# 538 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 110 "Parser.fsy"
                                      "?" 
                   )
# 110 "Parser.fsy"
                 : 'Modifier));
# 548 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 111 "Parser.fsy"
                                         "!" 
                   )
# 111 "Parser.fsy"
                 : 'Modifier));
# 558 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 112 "Parser.fsy"
                             "" 
                   )
# 112 "Parser.fsy"
                 : 'Modifier));
# 568 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 114 "Parser.fsy"
                                                    PRef(_1, None) 
                   )
# 114 "Parser.fsy"
                 : 'SimpleProduction));
# 579 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            let _3 = (let data = parseState.GetInput(3) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 115 "Parser.fsy"
                                                        
                             match ((_1), (_3)) with
                             | (("0",_), ("9", _)) -> makeToken ("NUMBER", (0,0)) "'0'..'9'"
                             | (("\\0",_), ("\\255", _)) -> makeToken ("CHAR", (0,0))  "'\\0'..'\\255'"
                             | (("a",_), ("z", _)) -> makeToken ("LOWER_LATIN", (0,0)) "'a'..'z'"
                             | (("A",_), ("Z", _)) -> makeToken ("UPPER_LATIN", (0,0)) "'A'..'Z'"
                             | ((a,_), (b, _)) -> makeToken ("", (0,0)) (a+".."+b)
                             
                   )
# 115 "Parser.fsy"
                 : 'SimpleProduction));
# 598 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 123 "Parser.fsy"
                                     PLiteral(_1) 
                   )
# 123 "Parser.fsy"
                 : 'SimpleProduction));
# 609 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 124 "Parser.fsy"
                                      PToken(_1) 
                   )
# 124 "Parser.fsy"
                 : 'SimpleProduction));
# 620 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'RulePart)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 127 "Parser.fsy"
                                             [_1] 
                   )
# 127 "Parser.fsy"
                 : 'RuleString));
# 631 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : 'RulePart)) in
            let _2 = (let data = parseState.GetInput(2) in (Microsoft.FSharp.Core.Operators.unbox data : 'RuleString)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 128 "Parser.fsy"
                                                  _1 :: _2 
                   )
# 128 "Parser.fsy"
                 : 'RuleString));
# 643 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 131 "Parser.fsy"
                                    ".." 
                   )
# 131 "Parser.fsy"
                 : 'RulePart));
# 653 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 132 "Parser.fsy"
                                 "~" 
                   )
# 132 "Parser.fsy"
                 : 'RulePart));
# 663 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 133 "Parser.fsy"
                                   "?" 
                   )
# 133 "Parser.fsy"
                 : 'RulePart));
# 673 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 134 "Parser.fsy"
                                       "!" 
                   )
# 134 "Parser.fsy"
                 : 'RulePart));
# 683 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 135 "Parser.fsy"
                                "+" 
                   )
# 135 "Parser.fsy"
                 : 'RulePart));
# 693 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 136 "Parser.fsy"
                                "*" 
                   )
# 136 "Parser.fsy"
                 : 'RulePart));
# 703 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 137 "Parser.fsy"
                                 "=" 
                   )
# 137 "Parser.fsy"
                 : 'RulePart));
# 713 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 138 "Parser.fsy"
                               "|" 
                   )
# 138 "Parser.fsy"
                 : 'RulePart));
# 723 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 139 "Parser.fsy"
                                  ")" 
                   )
# 139 "Parser.fsy"
                 : 'RulePart));
# 733 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 140 "Parser.fsy"
                                  "(" 
                   )
# 140 "Parser.fsy"
                 : 'RulePart));
# 743 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 141 "Parser.fsy"
                                   fst(_1) 
                   )
# 141 "Parser.fsy"
                 : 'RulePart));
# 754 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 142 "Parser.fsy"
                                      fst(_1) 
                   )
# 142 "Parser.fsy"
                 : 'RulePart));
# 765 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 143 "Parser.fsy"
                                    fst(_1) 
                   )
# 143 "Parser.fsy"
                 : 'RulePart));
# 776 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 144 "Parser.fsy"
                                       "{"+fst(_1)+"}" 
                   )
# 144 "Parser.fsy"
                 : 'RulePart));
# 787 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 145 "Parser.fsy"
                                       fst(_1) 
                   )
# 145 "Parser.fsy"
                 : 'RulePart));
# 798 "Parser.fs"
        (fun (parseState : Microsoft.FSharp.Text.Parsing.IParseState) ->
            let _1 = (let data = parseState.GetInput(1) in (Microsoft.FSharp.Core.Operators.unbox data : Source.t)) in
            Microsoft.FSharp.Core.Operators.box
                (
                   (
# 146 "Parser.fsy"
                                      fst(_1) 
                   )
# 146 "Parser.fsy"
                 : 'RulePart));
|]
# 810 "Parser.fs"
let tables () : Microsoft.FSharp.Text.Parsing.Tables<_> = 
  { reductions= _fsyacc_reductions ();
    endOfInputTag = _fsyacc_endOfInputTag;
    tagOfToken = tagOfToken;
    dataOfToken = _fsyacc_dataOfToken; 
    actionTableElements = _fsyacc_actionTableElements;
    actionTableRowOffsets = _fsyacc_actionTableRowOffsets;
    stateToProdIdxsTableElements = _fsyacc_stateToProdIdxsTableElements;
    stateToProdIdxsTableRowOffsets = _fsyacc_stateToProdIdxsTableRowOffsets;
    reductionSymbolCounts = _fsyacc_reductionSymbolCounts;
    immediateActions = _fsyacc_immediateActions;
    gotos = _fsyacc_gotos;
    sparseGotoTableRowOffsets = _fsyacc_sparseGotoTableRowOffsets;
    tagOfErrorTerminal = _fsyacc_tagOfErrorTerminal;
    parseError = (fun (ctxt:Microsoft.FSharp.Text.Parsing.ParseErrorContext<_>) -> 
                              match parse_error_rich with 
                              | Some f -> f ctxt
                              | None -> parse_error ctxt.Message);
    numTerminals = 26;
    productionToNonTerminalTable = _fsyacc_productionToNonTerminalTable  }
let engine lexer lexbuf startState = (tables ()).Interpret(lexer, lexbuf, startState)
let ParseAntlr lexer lexbuf : (Source.t, Source.t)Grammar.t * (string, string)System.Collections.Generic.Dictionary =
    Microsoft.FSharp.Core.Operators.unbox ((tables ()).Interpret(lexer, lexbuf, 0))

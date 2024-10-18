grammar AsymmetricCA;

program
    : statement* EOF
    ;

statement
    : caDeclaration ';'
    | cellDeclaration ';'
    | connection ';'
    | parameterSetting ';'
    | timeUnitSetting ';'
    | advanceSimulation
    | ruleDefinition
    ;

caDeclaration
    : 'CA' ID dimensions
    ;

dimensions
    : '[' INT (',' INT)* ']'
    ;

cellDeclaration
    : 'cell' ID cellPosition '{' cellProperties '}'
    ;

cellPosition
    : '[' INT (',' INT)* ']'
    ;

cellProperties
    : (propertyAssignment ';')*
    ;

propertyAssignment
    : ID '=' value
    ;

connection
    : 'connect' ID cellPosition 'with' ID cellPosition
    ;

parameterSetting
    : 'set' ID '=' value
    ;

timeUnitSetting
    : 'set' 'time_unit' '=' value
    ;

advanceSimulation
    : 'advance' INT 'steps' ';'
    | 'advance' 'step' ';'
    ;

ruleDefinition
    : 'rule' '{' ruleBody '}'
    ;

ruleBody
    : ( . )*?    // Match any character, non-greedy, until '}'
    ;

value
    : INT
    | FLOAT
    ;

ID  : [a-zA-Z_][a-zA-Z_0-9]*;
INT : [0-9]+;
FLOAT : [0-9]+'.'[0-9]+;

WS  : [ \t\r\n]+ -> skip;

// Skip comments outside the rule block
COMMENT : '#' ~[\r\n]* -> skip;

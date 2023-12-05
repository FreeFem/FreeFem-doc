// CodeMirror, copyright (c) by Marijn Haverbeke and others
// Distributed under an MIT license: https://codemirror.net/LICENSE

(function(mod) {
  if (typeof exports == "object" && typeof module == "object") // CommonJS
    mod(require("../../lib/codemirror"));
  else if (typeof define == "function" && define.amd) // AMD
    define(["../../lib/codemirror"], mod);
  else // Plain browser env
    mod(CodeMirror);
})(function(CodeMirror) {
"use strict";

function Context(indented, column, type, info, align, prev) {
  this.indented = indented;
  this.column = column;
  this.type = type;
  this.info = info;
  this.align = align;
  this.prev = prev;
}
function pushContext(state, col, type, info) {
  var indent = state.indented;
  if (state.context && state.context.type == "statement" && type != "statement")
    indent = state.context.indented;
  return state.context = new Context(indent, col, type, info, null, state.context);
}
function popContext(state) {
  var t = state.context.type;
  if (t == ")" || t == "]" || t == "}")
    state.indented = state.context.indented;
  return state.context = state.context.prev;
}

function typeBefore(stream, state, pos) {
  if (state.prevToken == "variable" || state.prevToken == "type") return true;
  if (/\S(?:[^- ]>|[*\]])\s*$|\*$/.test(stream.string.slice(0, pos))) return true;
  if (state.typeAtEndOfLine && stream.column() == stream.indentation()) return true;
}

function isTopScope(context) {
  for (;;) {
    if (!context || context.type == "top") return true;
    if (context.type == "}" && context.prev.info != "namespace") return false;
    context = context.prev;
  }
}

CodeMirror.defineMode("freefem", function(config, parserConfig) {
  var indentUnit = config.indentUnit,
      statementIndentUnit = parserConfig.statementIndentUnit || indentUnit,
      dontAlignCalls = parserConfig.dontAlignCalls,
      keywords = parserConfig.keywords || {},
      types = parserConfig.types || {},
      builtin = parserConfig.builtin || {},
      blockKeywords = parserConfig.blockKeywords || {},
      defKeywords = parserConfig.defKeywords || {},
      atoms = parserConfig.atoms || {},
      hooks = parserConfig.hooks || {},
      multiLineStrings = parserConfig.multiLineStrings,
      indentStatements = parserConfig.indentStatements !== false,
      indentSwitch = parserConfig.indentSwitch !== false,
      namespaceSeparator = parserConfig.namespaceSeparator,
      isPunctuationChar = parserConfig.isPunctuationChar || /[\[\]{}\(\),;\:\.]/,
      numberStart = parserConfig.numberStart || /[\d\.]/,
      number = parserConfig.number || /^(?:0x[a-f\d]+|0b[01]+|(?:\d+\.?\d*|\.\d+)(?:e[-+]?\d+)?)(u|ll?|l|f)?/i,
      isOperatorChar = parserConfig.isOperatorChar || /[#+'\-*&%=<>!?|\/]/,
      isIdentifierChar = parserConfig.isIdentifierChar || /[\w\$_\xa1-\uffff]/,
      // An optional function that takes a {string} token and returns true if it
      // should be treated as a builtin.
      isReservedIdentifier = parserConfig.isReservedIdentifier || false;

  var curPunc, isDefKeyword;

  function tokenBase(stream, state) {
    var ch = stream.next();
    if (hooks[ch]) {
      var result = hooks[ch](stream, state);
      if (result !== false) return result;
    }
    if (ch == '"' || ch == "'") {
      state.tokenize = tokenString(ch);
      return state.tokenize(stream, state);
    }
    if (isPunctuationChar.test(ch)) {
      curPunc = ch;
      return null;
    }
    if (numberStart.test(ch)) {
      stream.backUp(1)
      if (stream.match(number)) return "number"
      stream.next()
    }
    if (ch == "/") {
      if (stream.eat("*")) {
        state.tokenize = tokenComment;
        return tokenComment(stream, state);
      }
      if (stream.eat("/")) {
        stream.skipToEnd();
        return "comment";
      }
    }
    if (isOperatorChar.test(ch)) {
      while (!stream.match(/^\/[\/*]/, false) && stream.eat(isOperatorChar)) {}
      return "operator";
    }
    stream.eatWhile(isIdentifierChar);
    if (namespaceSeparator) while (stream.match(namespaceSeparator))
      stream.eatWhile(isIdentifierChar);

    var cur = stream.current();
    if (contains(keywords, cur)) {
      if (contains(blockKeywords, cur)) curPunc = "newstatement";
      if (contains(defKeywords, cur)) isDefKeyword = true;
      return "keyword";
    }
    if (contains(types, cur)) return "type";
    if (contains(builtin, cur)
        || (isReservedIdentifier && isReservedIdentifier(cur))) {
      if (contains(blockKeywords, cur)) curPunc = "newstatement";
      return "builtin";
    }
    if (contains(atoms, cur)) return "atom";
    return "variable";
  }

  function tokenString(quote) {
    return function(stream, state) {
      var escaped = false, next, end = false;
      while ((next = stream.next()) != null) {
        if (next == quote && !escaped) {end = true; break;}
        escaped = !escaped && next == "\\";
      }
      if (end || !(escaped || multiLineStrings))
        state.tokenize = null;
      return "string";
    };
  }

  function tokenComment(stream, state) {
    var maybeEnd = false, ch;
    while (ch = stream.next()) {
      if (ch == "/" && maybeEnd) {
        state.tokenize = null;
        break;
      }
      maybeEnd = (ch == "*");
    }
    return "comment";
  }

  function maybeEOL(stream, state) {
    if (parserConfig.typeFirstDefinitions && stream.eol() && isTopScope(state.context))
      state.typeAtEndOfLine = typeBefore(stream, state, stream.pos)
  }

  // Interface

  return {
    startState: function(basecolumn) {
      return {
        tokenize: null,
        context: new Context((basecolumn || 0) - indentUnit, 0, "top", null, false),
        indented: 0,
        startOfLine: true,
        prevToken: null
      };
    },

    token: function(stream, state) {
      var ctx = state.context;
      if (stream.sol()) {
        if (ctx.align == null) ctx.align = false;
        state.indented = stream.indentation();
        state.startOfLine = true;
      }
      if (stream.eatSpace()) { maybeEOL(stream, state); return null; }
      curPunc = isDefKeyword = null;
      var style = (state.tokenize || tokenBase)(stream, state);
      if (style == "comment" || style == "meta") return style;
      if (ctx.align == null) ctx.align = true;

      if (curPunc == ";" || curPunc == ":" || (curPunc == "," && stream.match(/^\s*(?:\/\/.*)?$/, false)))
        while (state.context.type == "statement") popContext(state);
      else if (curPunc == "{") pushContext(state, stream.column(), "}");
      else if (curPunc == "[") pushContext(state, stream.column(), "]");
      else if (curPunc == "(") pushContext(state, stream.column(), ")");
      else if (curPunc == "}") {
        while (ctx.type == "statement") ctx = popContext(state);
        if (ctx.type == "}") ctx = popContext(state);
        while (ctx.type == "statement") ctx = popContext(state);
      }
      else if (curPunc == ctx.type) popContext(state);
      else if (indentStatements &&
               (((ctx.type == "}" || ctx.type == "top") && curPunc != ";") ||
                (ctx.type == "statement" && curPunc == "newstatement"))) {
        pushContext(state, stream.column(), "statement", stream.current());
      }

      if (style == "variable" &&
          ((state.prevToken == "def" ||
            (parserConfig.typeFirstDefinitions && typeBefore(stream, state, stream.start) &&
             isTopScope(state.context) && stream.match(/^\s*\(/, false)))))
        style = "def";

      if (hooks.token) {
        var result = hooks.token(stream, state, style);
        if (result !== undefined) style = result;
      }

      if (style == "def" && parserConfig.styleDefs === false) style = "variable";

      state.startOfLine = false;
      state.prevToken = isDefKeyword ? "def" : style || curPunc;
      maybeEOL(stream, state);
      return style;
    },

    indent: function(state, textAfter) {
      if (state.tokenize != tokenBase && state.tokenize != null || state.typeAtEndOfLine) return CodeMirror.Pass;
      var ctx = state.context, firstChar = textAfter && textAfter.charAt(0);
      var closing = firstChar == ctx.type;
      if (ctx.type == "statement" && firstChar == "}") ctx = ctx.prev;
      if (parserConfig.dontIndentStatements)
        while (ctx.type == "statement" && parserConfig.dontIndentStatements.test(ctx.info))
          ctx = ctx.prev
      if (hooks.indent) {
        var hook = hooks.indent(state, ctx, textAfter, indentUnit);
        if (typeof hook == "number") return hook
      }
      var switchBlock = ctx.prev && ctx.prev.info == "switch";
      if (parserConfig.allmanIndentation && /[{(]/.test(firstChar)) {
        while (ctx.type != "top" && ctx.type != "}") ctx = ctx.prev
        return ctx.indented
      }
      if (ctx.type == "statement")
        return ctx.indented + (firstChar == "{" ? 0 : statementIndentUnit);
      if (ctx.align && (!dontAlignCalls || ctx.type != ")"))
        return ctx.column + (closing ? 0 : 1);
      if (ctx.type == ")" && !closing)
        return ctx.indented + statementIndentUnit;

      return ctx.indented + (closing ? 0 : indentUnit) +
        (!closing && switchBlock && !/^(?:case|default)\b/.test(textAfter) ? indentUnit : 0);
    },

    electricInput: indentSwitch ? /^\s*(?:case .*?:|default:|\{\}?|\})$/ : /^\s*[{}]$/,
    blockCommentStart: "/*",
    blockCommentEnd: "*/",
    blockCommentContinue: " * ",
    lineComment: "//",
    fold: "brace"
  };
});

  function words(str) {
    var obj = {}, words = str.split(" ");
    for (var i = 0; i < words.length; ++i) obj[words[i]] = true;
    return obj;
  }
  function contains(words, word) {
    if (typeof words === "function") {
      return words(word);
    } else {
      return words.propertyIsEnumerable(word);
    }
  }
  var ffKeywords = "break catch continue " +
          "else for if return try while";

  // Do not use this. Use the cTypes function below. This is global just to avoid
  // excessive calls when cTypes is being called multiple times during a parse.
  var ffTypes = "bool border complex dmatrix " +
          "fespace func gslspline ifstream " +
          "int macro matrix mesh mesh3 mpiComm " +
          "mpiGroup mpiRequest NewMacro EndMacro " +
          "ofstream Pmmap problem Psemaphore " +
          "real solve string varf";

  // Do not use this. Use the objCTypes function below. This is global just to avoid
  // excessive calls when objCTypes is being called multiple times during a parse.
  var basicObjCTypes = words("");

  // Returns true if identifier is a "C" type.
  // C type is defined as those that are reserved by the compiler (basicTypes),
  // and those that end in _t (Reserved by POSIX for types)
  // http://www.gnu.org/software/libc/manual/html_node/Reserved-Names.html
  function cTypes(identifier) {
    return contains(basicCTypes, identifier) || /.+_t/.test(identifier);
  }

  // Returns true if identifier is a "Objective C" type.
  function objCTypes(identifier) {
    return cTypes(identifier) || contains(basicObjCTypes, identifier);
  }

  var ffBlockKeywords = "catch else for if try while";
  var ffDefKeywords = "A A1 abserror absolute aniso aspectratio " +
          "B B1 bb beginend bin boundary bw " +
          "close cmm coef composante cutoff  " +
          "datafilename dataname dim distmax displacement doptions dparams " +
          "eps err errg " +
          "facemerge facetcl factorize file fill fixedborder flabel flags floatmesh floatsol fregion " +
          "gradation grey " +
          "hmax hmin holelist hsv  " +
          "init inquire inside IsMetric iso ivalue " +
          "keepbackvertices " +
          "label labeldown labelmid labelup levelset loptions lparams " +
          "maxit maxsubdiv meditff mem memory metric mode " +
          "nbarrow nbiso nbiter nbjacoby nboffacetcl nbofholes nbofregions nbregul nbsmooth nbvx ncv nev nomeshgeneration normalization " +
          "omega op optimize option options order orientation " +
          "periodic power precon prev ps ptmerge " +
          "qfe qforder qft qfV " +
          "ratio rawvector reffacelow reffacemid reffaceup refnum reftet reftri region regionlist renumv rescaling ridgeangle " +
          "save sigma sizeofvolume smoothing solver sparams split splitin2 splitpbedge stop strategy swap switch sym " +
          "t tgv thetamax tol tolpivot tolpivotsym transfo " +
          "U2Vc " +
          "value varrow vector veps viso " +
          "wait width withsurfacemesh WindowIndex which " +
          "zbound";

  function cppHook(stream, state) {
    if (!state.startOfLine) return false
    for (var ch, next = null; ch = stream.peek();) {
      if (ch == "\\" && stream.match(/^.$/)) {
        next = cppHook
        break
      } else if (ch == "/" && stream.match(/^\/[\/\*]/, false)) {
        break
      }
      stream.next()
    }
    state.tokenize = next
    return "meta"
  }

  function pointerHook(_stream, state) {
    if (state.prevToken == "type") return "type";
    return false;
  }

  function ffIsReservedIdentifier(token) {
    if (!token || token.length < 2) return false;
	var list = words("BDM1 BDM1Ortho Edge03d " +
            "Edge13d Edge23d FEQF HCT " +
            "P0 P03d P0Edge P1 P13d P1b " +
            "P1b3d P1bl P1bl3d P1dc P1Edge " +
            "P1nc P2 P23d P2b P2BR P2dc " +
            "P2Edge P2h P2Morley P2pnc " +
            "P3 P3dc P3Edge P4 P4dc P4Edge " +
            "P5Edge RT0 RT03d RT0Ortho RT1 " +
            "RT1Ortho RT2 RT2Ortho " +

            "qf1pE qf1pElump qf1pT qf1pTlump " +
            "qfV1 qfV1lump qf2pE qf2pT qf2pT4P1 " +
            "qfV2 qf3pE qf4pE qf5pE qf5pT qfV5 " +
            "qf7pT qf9pT qfnbpE " +

            "ARGV append area be binary BoundaryEdge bordermeasure  " +
            "CG Cholesky cin cout Crout default diag edgeOrientation  " +
            "endl FILE fixed GMRES good hTriangle im imax imin InternalEdge  " +
            "l1 l2 label lenEdge length LINE linfty LU m max measure min " +
            "mpiAnySource mpiBAND mpiBXOR mpiCommWorld mpiLAND mpiLOR " +
            "mpiLXOR mpiMAX mpiMIN mpiPROD mpirank mpisize mpiSUM mpiUndefined " +
            "n N nbe ndof ndofK noshowbase noshowpos notaregion nt nTonEdge " +
            "nuEdge nuTriangle nv P pi precision quantile re region " +
            "scientific searchMethod setw showbase showpos sparsesolver " +
            "sum tellp true UMFPACK unused whoinElement verbosity version " +
            "volume x y z");
	if (Object.keys(list).includes(token))
		return true
	return false
  }

  var ffBuiltIn = "abs acos acosh adaptmesh adj AffineCG AffineGMRES arg asin asinh assert atan atan2 atanh atof atoi " +
          "BFGS broadcast buildlayers buildmesh " +
          "ceil chi complexEigenValue copysign change checkmovemesh clock cmaes conj convect cos cosh cube " +
          "d dd dfft diffnp diffpos dimKrylov dist dumptable dx dxx dxy dxz dy dyx dyy dyz dz dzx dzy dzz " +
          "EigenValue emptymesh erf erfc exec exit exp " +
          "fdim ffind find floor flush fmax fmin fmod freeyams " +
          "getARGV getline gmshload gmshload3 " +
          "gslcdfugaussianP gslcdfugaussianQ gslcdfugaussianPinv gslcdfugaussianQinv gslcdfgaussianP gslcdfgaussianQ " +
          "gslcdfgaussianPinv gslcdfgaussianQinv gslcdfgammaP gslcdfgammaQ gslcdfgammaPinv gslcdfgammaQinv gslcdfcauchyP " +
          "gslcdfcauchyQ gslcdfcauchyPinv gslcdfcauchyQinv gslcdflaplaceP gslcdflaplaceQ gslcdflaplacePinv gslcdflaplaceQinv " +
          "gslcdfrayleighP gslcdfrayleighQ gslcdfrayleighPinv gslcdfrayleighQinv gslcdfchisqP gslcdfchisqQ gslcdfchisqPinv " +
          "gslcdfchisqQinv gslcdfexponentialP gslcdfexponentialQ gslcdfexponentialPinv gslcdfexponentialQinv gslcdfexppowP " +
          "gslcdfexppowQ gslcdftdistP gslcdftdistQ gslcdftdistPinv gslcdftdistQinv gslcdffdistP gslcdffdistQ gslcdffdistPinv " +
          "gslcdffdistQinv gslcdfbetaP gslcdfbetaQ gslcdfbetaPinv gslcdfbetaQinv gslcdfflatP gslcdfflatQ gslcdfflatPinv gslcdfflatQinv " +
          "gslcdflognormalP gslcdflognormalQ gslcdflognormalPinv gslcdflognormalQinv gslcdfgumbel1P gslcdfgumbel1Q gslcdfgumbel1Pinv " +
          "gslcdfgumbel1Qinv gslcdfgumbel2P gslcdfgumbel2Q gslcdfgumbel2Pinv gslcdfgumbel2Qinv gslcdfweibullP gslcdfweibullQ " +
          "gslcdfweibullPinv gslcdfweibullQinv gslcdfparetoP gslcdfparetoQ gslcdfparetoPinv gslcdfparetoQinv gslcdflogisticP " +
          "gslcdflogisticQ gslcdflogisticPinv gslcdflogisticQinv gslcdfbinomialP gslcdfbinomialQ gslcdfpoissonP gslcdfpoissonQ " +
          "gslcdfgeometricP gslcdfgeometricQ gslcdfnegativebinomialP gslcdfnegativebinomialQ gslcdfpascalP gslcdfpascalQ " +
          "gslinterpakima gslinterpakimaperiodic gslinterpcsplineperiodic gslinterpcspline gslinterpsteffen gslinterplinear " +
          "gslinterppolynomial gslranbernoullipdf gslranbeta gslranbetapdf gslranbinomialpdf gslranexponential gslranexponentialpdf " +
          "gslranexppow gslranexppowpdf gslrancauchy gslrancauchypdf gslranchisq gslranchisqpdf gslranerlang gslranerlangpdf " +
          "gslranfdist gslranfdistpdf gslranflat gslranflatpdf gslrangamma gslrangammaint gslrangammapdf gslrangammamt " +
          "gslrangammaknuth gslrangaussian gslrangaussianratiomethod gslrangaussianziggurat gslrangaussianpdf gslranugaussian " +
          "gslranugaussianratiomethod gslranugaussianpdf gslrangaussiantail gslrangaussiantailpdf gslranugaussiantail " +
          "gslranugaussiantailpdf gslranlandau gslranlandaupdf gslrangeometricpdf gslrangumbel1 gslrangumbel1pdf gslrangumbel2 " +
          "gslrangumbel2pdf gslranlogistic gslranlogisticpdf gslranlognormal gslranlognormalpdf gslranlogarithmicpdf " +
          "gslrannegativebinomialpdf gslranpascalpdf gslranpareto gslranparetopdf gslranpoissonpdf gslranrayleigh " +
          "gslranrayleighpdf gslranrayleightail gslranrayleightailpdf gslrantdist gslrantdistpdf gslranlaplace gslranlaplacepdf " +
          "gslranlevy gslranweibull gslranweibullpdf gslsfairyAi gslsfairyBi gslsfairyAiscaled gslsfairyBiscaled " +
          "gslsfairyAideriv gslsfairyBideriv gslsfairyAiderivscaled gslsfairyBiderivscaled gslsfairyzeroAi gslsfairyzeroBi " +
          "gslsfairyzeroAideriv gslsfairyzeroBideriv gslsfbesselJ0 gslsfbesselJ1 gslsfbesselJn gslsfbesselY0 gslsfbesselY1 " +
          "gslsfbesselYn gslsfbesselI0 gslsfbesselI1 gslsfbesselIn gslsfbesselI0scaled gslsfbesselI1scaled gslsfbesselInscaled " +
          "gslsfbesselK0 gslsfbesselK1 gslsfbesselKn gslsfbesselK0scaled gslsfbesselK1scaled gslsfbesselKnscaled gslsfbesselj0 " +
          "gslsfbesselj1 gslsfbesselj2 gslsfbesseljl gslsfbessely0 gslsfbessely1 gslsfbessely2 gslsfbesselyl gslsfbesseli0scaled " +
          "gslsfbesseli1scaled gslsfbesseli2scaled gslsfbesselilscaled gslsfbesselk0scaled gslsfbesselk1scaled gslsfbesselk2scaled " +
          "gslsfbesselklscaled gslsfbesselJnu gslsfbesselYnu gslsfbesselInuscaled gslsfbesselInu gslsfbesselKnuscaled gslsfbesselKnu " +
          "gslsfbessellnKnu gslsfbesselzeroJ0 gslsfbesselzeroJ1 gslsfbesselzeroJnu gslsfclausen gslsfhydrogenicR1 gslsfdawson " +
          "gslsfdebye1 gslsfdebye2 gslsfdebye3 gslsfdebye4 gslsfdebye5 gslsfdebye6 gslsfdilog gslsfmultiply gslsfellintKcomp " +
          "gslsfellintEcomp gslsfellintPcomp gslsfellintDcomp gslsfellintF gslsfellintE gslsfellintRC gslsferfc gslsflogerfc gslsferf " +
          "gslsferfZ gslsferfQ gslsfhazard gslsfexp gslsfexpmult gslsfexpm1 gslsfexprel gslsfexprel2 gslsfexpreln gslsfexpintE1 " +
          "gslsfexpintE2 gslsfexpintEn gslsfexpintE1scaled gslsfexpintE2scaled gslsfexpintEnscaled gslsfexpintEi gslsfexpintEiscaled " +
          "gslsfShi gslsfChi gslsfexpint3 gslsfSi gslsfCi gslsfatanint gslsffermidiracm1 gslsffermidirac0 gslsffermidirac1 " +
          "gslsffermidirac2 gslsffermidiracint gslsffermidiracmhalf gslsffermidirachalf gslsffermidirac3half gslsffermidiracinc0 " +
          "gslsflngamma gslsfgamma gslsfgammastar gslsfgammainv gslsftaylorcoeff gslsffact gslsfdoublefact gslsflnfact gslsflndoublefact " +
          "gslsflnchoose gslsfchoose gslsflnpoch gslsfpoch gslsfpochrel gslsfgammaincQ gslsfgammaincP gslsfgammainc gslsflnbeta " +
          "gslsfbeta gslsfbetainc gslsfgegenpoly1 gslsfgegenpoly2 gslsfgegenpoly3 gslsfgegenpolyn gslsfhyperg0F1 gslsfhyperg1F1int " +
          "gslsfhyperg1F1 gslsfhypergUint gslsfhypergU gslsfhyperg2F0 gslsflaguerre1 gslsflaguerre2 gslsflaguerre3 gslsflaguerren " +
          "gslsflambertW0 gslsflambertWm1 gslsflegendrePl gslsflegendreP1 gslsflegendreP2 gslsflegendreP3 gslsflegendreQ0 gslsflegendreQ1 " +
          "gslsflegendreQl gslsflegendrePlm gslsflegendresphPlm gslsflegendrearraysize gslsfconicalPhalf gslsfconicalPmhalf gslsfconicalP0 " +
          "gslsfconicalP1 gslsfconicalPsphreg gslsfconicalPcylreg gslsflegendreH3d0 gslsflegendreH3d1 gslsflegendreH3d gslsflog gslsflogabs " +
          "gslsflog1plusx gslsflog1plusxmx gslsfpowint gslsfpsiint gslsfpsi gslsfpsi1piy gslsfpsi1int gslsfpsi1 gslsfpsin " +
          "gslsfsynchrotron1 gslsfsynchrotron2 gslsftransport2 gslsftransport3 gslsftransport4 gslsftransport5 gslsfsin gslsfcos " +
          "gslsfhypot gslsfsinc gslsflnsinh gslsflncosh gslsfanglerestrictsymm gslsfanglerestrictpos gslsfzetaint gslsfzeta gslsfzetam1 " +
          "gslsfzetam1int gslsfhzeta gslsfetaint gslsfeta " +
          "imag int1d int2d int3d intalledges intallfaces interpolate invdiff invdiffnp invdiffpos Isend isInf isNaN isoline Irecv " +
          "j0 j1 jn jump " +
          "lgamma LinearCG LinearGMRES log log10 lrint lround " +
          "max mean medit min mmg3d movemesh movemesh23 mpiAlltoall mpiAlltoallv mpiAllgather mpiAllgatherv mpiAllReduce mpiBarrier " +
          "mpiGather mpiGatherv mpiRank mpiReduce mpiScatter mpiScatterv mpiSize mpiWait mpiWaitAny mpiWtick mpiWtime mshmet " +
          "NLCG " +
          "on " +
          "plot polar Post pow processor processorblock projection " +
          "randinit randint31 randint32 random randreal1 randreal2 randreal3 randres53 Read readmesh readmesh3 Recv rfind rint round " +
          "savemesh savesol savevtk seekg Sent set sign signbit sin sinh sort splitComm splitmesh sqrt square srandom srandomdev Stringification swap system " +
          "tan tanh tellg tetg tetgconvexhull tetgreconstruction tetgtransfo tgamma triangulate trunc " +
          "Wait Write " +
          "y0 y1 yn";

  function cpp14Literal(stream) {
    stream.eatWhile(/[\w\.']/);
    return "number";
  }

  function cpp11StringHook(stream, state) {
    stream.backUp(1);
    // Raw strings.
    if (stream.match(/(R|u8R|uR|UR|LR)/)) {
      var match = stream.match(/"([^\s\\()]{0,16})\(/);
      if (!match) {
        return false;
      }
      state.cpp11RawStringDelim = match[1];
      state.tokenize = tokenRawString;
      return tokenRawString(stream, state);
    }
    // Unicode strings/chars.
    if (stream.match(/(u8|u|U|L)/)) {
      if (stream.match(/["']/, /* eat */ false)) {
        return "string";
      }
      return false;
    }
    // Ignore this hook.
    stream.next();
    return false;
  }

  function cppLooksLikeConstructor(word) {
    var lastTwo = /(\w+)::~?(\w+)$/.exec(word);
    return lastTwo && lastTwo[1] == lastTwo[2];
  }

  // C#-style strings where "" escapes a quote.
  function tokenAtString(stream, state) {
    var next;
    while ((next = stream.next()) != null) {
      if (next == '"' && !stream.eat('"')) {
        state.tokenize = null;
        break;
      }
    }
    return "string";
  }

  // C++11 raw string literal is <prefix>"<delim>( anything )<delim>", where
  // <delim> can be a string up to 16 characters long.
  function tokenRawString(stream, state) {
    // Escape characters that have special regex meanings.
    var delim = state.cpp11RawStringDelim.replace(/[^\w\s]/g, '\\$&');
    var match = stream.match(new RegExp(".*?\\)" + delim + '"'));
    if (match)
      state.tokenize = null;
    else
      stream.skipToEnd();
    return "string";
  }

  function def(mimes, mode) {
    if (typeof mimes == "string") mimes = [mimes];
    var words = [];
    function add(obj) {
      if (obj) for (var prop in obj) if (obj.hasOwnProperty(prop))
        words.push(prop);
    }
    add(mode.keywords);
    add(mode.types);
    add(mode.builtin);
    add(mode.atoms);
    if (words.length) {
      mode.helperType = mimes[0];
      CodeMirror.registerHelper("hintWords", mimes[0], words);
    }

    for (var i = 0; i < mimes.length; ++i)
      CodeMirror.defineMIME(mimes[i], mode);
  }

  def(["text/x-ff++src"], {
    name: "freefem",
    // Keywords from https://en.cppreference.com/w/cpp/keyword includes C++20.
    keywords: words(ffKeywords),
    types: words(ffTypes),
    blockKeywords: words(ffBlockKeywords),
	builtin: words(ffBuiltIn),
	defKeywords: words(ffDefKeywords),
    typeFirstDefinitions: true,
    dontIndentStatements: /^template$/,
    isIdentifierChar: /[\w\$_~\xa1-\uffff]/,
    isReservedIdentifier: ffIsReservedIdentifier,
    hooks: {
      "0": cpp14Literal,
      "1": cpp14Literal,
      "2": cpp14Literal,
      "3": cpp14Literal,
      "4": cpp14Literal,
      "5": cpp14Literal,
      "6": cpp14Literal,
      "7": cpp14Literal,
      "8": cpp14Literal,
      "9": cpp14Literal,
      token: function(stream, state, style) {
        if (style == "variable" && stream.peek() == "(" &&
            (state.prevToken == ";" || state.prevToken == null ||
             state.prevToken == "}") &&
            cppLooksLikeConstructor(stream.current()))
          return "def";
      }
    },
    namespaceSeparator: "::",
    modeProps: {fold: ["brace", "include"]}
  });
});

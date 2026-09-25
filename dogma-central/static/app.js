"use strict";

// ---------------------------------------------------------------------------
// Datos de referencia
// ---------------------------------------------------------------------------
const ENZIMAS_REPLICACION = {
  iniciadoras: { corto: "Iniciadoras", nombre: "Proteínas iniciadoras", funcion: "Reconocen el origen de replicación y cargan la helicasa." },
  helicasa: { corto: "Hel", nombre: "Helicasa", funcion: "Rompe los puentes de hidrógeno entre bases y abre la doble hélice. Gasta ATP." },
  topoisomerasa: { corto: "Topo", nombre: "Topoisomerasa (girasa)", funcion: "Corta y vuelve a unir el ADN por delante de la horquilla para aliviar el superenrollamiento." },
  ssb: { corto: "SSB", nombre: "Proteínas SSB", funcion: "Se unen al ADN de cadena sencilla e impiden que se vuelva a aparear." },
  primasa: { corto: "Primasa", nombre: "Primasa", funcion: "ARN polimerasa que fabrica los cebadores de ARN sobre los que trabaja la ADN polimerasa." },
  pol3: { corto: "Pol III", nombre: "ADN polimerasa III", funcion: "Alarga las cadenas nuevas en sentido 5'→3' desde el 3'-OH de un cebador y corrige sus propios errores." },
  pol1: { corto: "Pol I", nombre: "ADN polimerasa I", funcion: "Retira los cebadores de ARN con su actividad exonucleasa 5'→3' y los sustituye por ADN." },
  ligasa: { corto: "Ligasa", nombre: "ADN ligasa", funcion: "Sella las mellas con un enlace fosfodiéster y une los fragmentos de Okazaki." },
};

const ETAPAS = ["replicacion", "transcripcion", "traduccion"];

const estado = {
  datos: null,
  etapa: "replicacion",
  paso: { replicacion: 0, transcripcion: 0, traduccion: 0 },
  temporizador: null,
};

const $ = (id) => document.getElementById(id);

function el(etiqueta, clase, texto) {
  const e = document.createElement(etiqueta);
  if (clase) e.className = clase;
  if (texto !== undefined && texto !== null) e.textContent = texto;
  return e;
}

function secuenciaColoreada(seq, claseExtra) {
  const cont = el("span", "seq " + (claseExtra || ""));
  for (const b of seq) cont.appendChild(el("span", b, b));
  return cont;
}

// ---------------------------------------------------------------------------
// Rejilla de secuencias (una columna por nucleótido)
// ---------------------------------------------------------------------------
function crearFila(n, nombre, izq, der, clase) {
  const f = el("div", "fila " + (clase || ""));
  f.style.setProperty("--n", n);
  const nom = el("div", "nombre");
  nom.style.gridColumn = "1";
  nom.appendChild(el("span", "", nombre));
  if (izq) nom.appendChild(el("span", "extremo-izq", izq));
  f.appendChild(nom);
  const fin = el("div", "extremo", der || "");
  fin.style.gridColumn = String(n + 2);
  f.appendChild(fin);
  return f;
}

function celda(fila, i, texto, clases, titulo) {
  const c = el("div", "b " + (clases || ""), texto);
  c.style.gridColumn = String(i + 2);
  if (titulo) c.title = titulo;
  fila.appendChild(c);
  return c;
}

function chip(fila, clave, inicio, fin, etiqueta, titulo) {
  const c = el("div", "chip " + clave, etiqueta);
  c.style.gridColumn = `${inicio + 2} / ${Math.max(fin, inicio + 1) + 2}`;
  c.title = titulo || etiqueta;
  fila.appendChild(c);
}

// ---------------------------------------------------------------------------
// Replicación
// ---------------------------------------------------------------------------
function dibujarReplicacion(rep, foto) {
  const n = rep.superior.length;
  const h = foto.horquilla;
  const mellas = new Set(foto.mellas);
  const rej = el("div", "rejilla");

  const filaEnzimas = (hebra) => {
    const f = crearFila(n, "", "", "", "enzimas");
    for (const e of foto.enzimas.filter((x) => x.hebra === hebra)) {
      const info = ENZIMAS_REPLICACION[e.clave];
      chip(f, e.clave, e.inicio, e.fin, info.corto, info.nombre);
    }
    return f;
  };

  const filaNueva = (hebra, nombre, izq, der) => {
    const datos = hebra === "lider" ? foto.lider : foto.rezagada;
    const f = crearFila(n, nombre, izq, der, "nueva");
    for (let i = 0; i < n; i++) {
      if (i >= h) { celda(f, i, "", "enlace"); continue; }
      const tipo = datos.tipos[i];
      const base = datos.bases[i];
      if (tipo === "-") { celda(f, i, "·", "vacia"); continue; }
      let clases = "b-" + base;
      let titulo;
      if (tipo === "r") { clases += " cebador"; titulo = "Cebador de ARN"; }
      else if (hebra === "lider") { clases += " lider"; titulo = "Cadena líder (ADN)"; }
      else {
        const fr = datos.fragmentos[i];
        clases += fr % 2 ? " frag-1" : " frag-2";
        titulo = `Fragmento de Okazaki ${fr}`;
      }
      if (hebra === "rezagada") {
        if (mellas.has(i)) clases += " mella";
        const fr = datos.fragmentos[i];
        if (i > 0 && datos.tipos[i - 1] !== "-" && datos.fragmentos[i - 1] !== fr) clases += " inicio-frag";
      }
      celda(f, i, base, clases, titulo);
    }
    return f;
  };

  const filaParental = (seq, nombre, izq, der, nuevaTipos) => {
    const f = crearFila(n, nombre, izq, der, "parental");
    for (let i = 0; i < n; i++) {
      const expuesta = i < h && nuevaTipos[i] === "-";
      celda(f, i, seq[i], expuesta ? "expuesta" : "", expuesta ? "Cadena sencilla expuesta" : "");
    }
    return f;
  };

  const eje = crearFila(n, "", "", "", "eje");
  for (let i = 0; i < n; i++) celda(eje, i, "", i >= h ? "enlace" : "");
  for (const e of foto.enzimas.filter((x) => x.hebra === "horquilla")) {
    const info = ENZIMAS_REPLICACION[e.clave];
    chip(eje, e.clave, e.inicio, e.fin, info.corto, info.nombre);
  }

  rej.append(
    filaEnzimas("rezagada"),
    filaParental(rep.superior, "Parental superior", "5'", "3'", foto.rezagada.tipos),
    filaNueva("rezagada", "Nueva: rezagada", "3'", "5'"),
    eje,
    filaNueva("lider", "Nueva: líder", "5'", "3'"),
    filaParental(rep.inferior, "Parental inferior", "3'", "5'", foto.lider.tipos),
    filaEnzimas("lider"),
  );
  return rej;
}

function leyendaReplicacion(foto) {
  const cont = el("div");
  const muestras = el("div", "leyenda");
  const items = [
    ["cebador", "Cebador de ARN"],
    ["frag", "Fragmentos de Okazaki (alternos)"],
    ["lider", "Cadena líder"],
    ["mella", "Mella"],
    ["enlace", "Puentes de hidrógeno"],
  ];
  for (const [clave, texto] of items) {
    const it = el("span", "item");
    const m = el("span", "muestra " + clave);
    if (clave === "frag") m.style.background = "linear-gradient(90deg, var(--frag-1) 50%, var(--frag-2) 50%)";
    if (clave === "lider") m.style.background = "var(--lider)";
    it.append(m, document.createTextNode(texto));
    muestras.appendChild(it);
  }
  cont.appendChild(muestras);

  const activas = new Set(foto.enzimas.map((e) => e.clave));
  const lista = el("div", "enzimas-lista");
  for (const [clave, info] of Object.entries(ENZIMAS_REPLICACION)) {
    const d = el("div", activas.has(clave) ? "activa" : "");
    d.appendChild(el("strong", "", info.nombre + ". "));
    d.appendChild(document.createTextNode(info.funcion));
    lista.appendChild(d);
  }
  cont.appendChild(lista);
  cont.appendChild(el("p", "nota",
    "Nombres de las enzimas de E. coli. En eucariotas, el cebador lo pone el complejo Pol α-primasa, " +
    "alargan las polimerasas δ (rezagada) y ε (líder), y los cebadores los retiran la RNasa H y FEN1."));
  return cont;
}

function extraReplicacion(rep, foto) {
  const cont = el("div");
  const frags = rep.fragmentos
    .map((f) => `${f.numero} (${f.inicio + 1}–${f.fin})`)
    .join(", ");
  cont.appendChild(el("p", "nota", `Fragmentos de Okazaki, por posiciones: ${frags}.`));

  if (foto.etapa === "resultado") {
    const hijas = el("div", "hijas");
    rep.hijas.forEach((h) => {
      const d = el("div", "hija");
      d.appendChild(el("h4", "", `${h.nombre}: conserva la hebra parental ${h.parental}`));
      const linea = (etq, seq, nueva) => {
        const p = el("div", "seq");
        p.appendChild(document.createTextNode(etq + "-"));
        const s = secuenciaColoreada(seq, nueva ? "nueva" : "");
        p.appendChild(s);
        p.appendChild(document.createTextNode(etq === "5'" ? "-3'" : "-5'"));
        return p;
      };
      d.appendChild(linea("5'", h.superior, h.parental !== "superior"));
      d.appendChild(linea("3'", h.inferior, h.parental !== "inferior"));
      hijas.appendChild(d);
    });
    cont.appendChild(hijas);
    cont.appendChild(el("p", rep.copia_fiel ? "nota ok" : "nota aviso",
      rep.copia_fiel
        ? "Las dos moléculas hijas son idénticas a la parental. La hebra nueva de cada una va sombreada."
        : "Las moléculas hijas no coinciden con la parental."));
  }
  return cont;
}

// ---------------------------------------------------------------------------
// Transcripción
// ---------------------------------------------------------------------------
const HIBRIDO = 8; // longitud aproximada del híbrido ARN-ADN en la burbuja

function dibujarTranscripcion(tr, foto) {
  const n = tr.superior.length;
  const k = foto.transcritos;
  const derecha = tr.sentido === "derecha";
  const terminada = foto.etapa === "terminacion";

  const transcrito = (i) => (derecha ? i < k : i >= n - k);
  let bIni = 0, bFin = 0, pol = -1;
  if (!terminada) {
    if (derecha) { bIni = Math.max(0, k - HIBRIDO); bFin = Math.min(n, k + 3); pol = Math.min(k, n - 1); }
    else { bIni = Math.max(0, n - k - 3); bFin = Math.min(n, n - k + HIBRIDO); pol = Math.max(n - 1 - k, 0); }
  }
  const enBurbuja = (i) => i >= bIni && i < bFin;

  const rej = el("div", "rejilla");

  const filaADN = (seq, nombre, izq, der) => {
    const f = crearFila(n, nombre, izq, der, "parental");
    for (let i = 0; i < n; i++) celda(f, i, seq[i], enBurbuja(i) ? "expuesta" : "");
    return f;
  };
  const eje = crearFila(n, "", "", "", "eje");
  for (let i = 0; i < n; i++) celda(eje, i, "", enBurbuja(i) ? "burbuja" : "enlace");

  const filaARN = derecha
    ? crearFila(n, "ARNm", "5'", "3'", "nueva")
    : crearFila(n, "ARNm", "3'", "5'", "nueva");
  for (let i = 0; i < n; i++) {
    if (!transcrito(i)) { celda(filaARN, i, "", ""); continue; }
    const b = tr.arn_alineado[i];
    const clase = "b-" + b + (enBurbuja(i) ? " burbuja" : " arn-libre");
    celda(filaARN, i, b, clase, enBurbuja(i) ? "Híbrido ARN-ADN" : "ARN ya separado del molde");
  }

  const filaEnz = crearFila(n, "", "", "", "enzimas");
  if (pol >= 0) {
    const ini = Math.max(0, pol - 1);
    chip(filaEnz, "arnpol", ini, Math.min(n, ini + 3), "ARN pol", "ARN polimerasa");
  }

  const sup = filaADN(tr.superior, tr.hebra_molde === "superior" ? "Molde (superior)" : "Codificante (superior)", "5'", "3'");
  const inf = filaADN(tr.inferior, tr.hebra_molde === "inferior" ? "Molde (inferior)" : "Codificante (inferior)", "3'", "5'");

  if (derecha) rej.append(filaEnz, sup, eje, inf, filaARN);
  else rej.append(filaARN, sup, eje, inf, filaEnz);
  return rej;
}

function leyendaTranscripcion() {
  const cont = el("div");
  const muestras = el("div", "leyenda");
  const add = (estilo, texto, clase) => {
    const it = el("span", "item");
    const m = el("span", "muestra " + (clase || ""));
    if (estilo) m.style.background = estilo;
    it.append(m, document.createTextNode(texto));
    muestras.appendChild(it);
  };
  add("#fff4cc", "Burbuja de transcripción (híbrido ARN-ADN)");
  add("#e8ddef", "ARN ya separado del molde");
  add(null, "Puentes de hidrógeno", "enlace");
  cont.appendChild(muestras);
  const lista = el("div", "enzimas-lista");
  const d = el("div", "activa");
  d.appendChild(el("strong", "", "ARN polimerasa. "));
  d.appendChild(document.createTextNode(
    "Abre la doble hélice, lee la hebra molde 3'→5' y une ribonucleótidos en sentido 5'→3'. No necesita cebador."));
  lista.appendChild(d);
  const d2 = el("div", "");
  d2.appendChild(el("strong", "", "Reglas de emparejamiento. "));
  d2.appendChild(document.createTextNode("Molde A → U, T → A, C → G, G → C en el ARN."));
  lista.appendChild(d2);
  cont.appendChild(lista);
  return cont;
}

function extraTranscripcion(tr, foto) {
  const cont = el("div", "arn-lineal");
  cont.appendChild(el("span", "rot", "ARNm"));
  cont.appendChild(document.createTextNode("5'-"));
  cont.appendChild(secuenciaColoreada(tr.arnm.slice(0, foto.transcritos)));
  if (foto.transcritos < tr.arnm.length) cont.appendChild(el("span", "pendiente", "…"));
  cont.appendChild(document.createTextNode("-3'"));
  return cont;
}

// ---------------------------------------------------------------------------
// Traducción
// ---------------------------------------------------------------------------
function dibujarTraduccion(td, foto) {
  const cont = el("div");
  const pista = el("div", "pista");
  pista.appendChild(el("span", "extremo-arn", "5'"));

  if (td.inicio === -1) {
    pista.appendChild(el("span", "utr", td.arnm));
    pista.appendChild(el("span", "extremo-arn", "3'"));
    cont.appendChild(pista);
    return cont;
  }

  if (td.inicio > 0) pista.appendChild(el("span", "utr", td.arnm.slice(0, td.inicio)));

  const s = foto.sitios;
  const actual = s.A ?? s.P ?? 0;
  td.codones.forEach((c, j) => {
    const cod = el("div", "codon");
    if (j <= actual) cod.classList.add("leido");
    if (c.parada) cod.classList.add("es-parada");
    let sitio = null;
    for (const nombre of ["E", "P", "A"]) if (s[nombre] === j) sitio = nombre;
    if (sitio) {
      cod.classList.add("sitio", "sitio-" + sitio);
      if (sitio === "P" && s.E === null) cod.classList.add("primero");
      if (sitio === "P" && s.A === null) cod.classList.add("ultimo");
      cod.appendChild(el("span", "rotulo-sitio", sitio));
    }
    const tri = el("div", "tri");
    for (const b of c.codon) tri.appendChild(el("span", b, b));
    cod.appendChild(tri);
    cod.appendChild(el("span", "num-codon", String(j + 1)));

    const arnt = el("div", "arnt");
    if (sitio === "P") {
      arnt.appendChild(el("span", "anti", c.anticodon));
      arnt.appendChild(el("span", "aa-bola", c.abrev));
      arnt.title = `ARNt con anticodón 3'-${c.anticodon}-5' unido a la cadena en crecimiento`;
    } else if (sitio === "E") {
      arnt.classList.add("vacio");
      arnt.appendChild(el("span", "anti", c.anticodon));
      arnt.appendChild(el("span", "aa-bola", "sale"));
      arnt.title = "ARNt descargado que abandona el ribosoma";
    } else if (sitio === "A") {
      if (c.parada && foto.etapa === "terminacion") {
        arnt.classList.add("factor");
        arnt.appendChild(el("span", "anti", " "));
        arnt.appendChild(el("span", "aa-bola", "Factor"));
        arnt.title = "Factor de liberación";
      } else {
        arnt.classList.add("vacio");
        arnt.appendChild(el("span", "anti", " "));
        arnt.appendChild(el("span", "aa-bola", "libre"));
        arnt.title = "Sitio A libre, a la espera del siguiente aminoacil-ARNt";
      }
    }
    cod.appendChild(arnt);
    pista.appendChild(cod);
  });

  const finMarco = td.inicio + td.codones.length * 3;
  if (finMarco < td.arnm.length) pista.appendChild(el("span", "utr", td.arnm.slice(finMarco)));
  pista.appendChild(el("span", "extremo-arn", "3'"));
  cont.appendChild(pista);

  // Cadena polipeptídica
  const cadena = el("div", "cadena");
  const liberada = foto.etapa === "terminacion";
  cadena.appendChild(el("h4", "", liberada
    ? `Proteína liberada: ${foto.peptido.length} aminoácidos`
    : `Cadena en crecimiento: ${foto.peptido.length} aminoácido${foto.peptido.length === 1 ? "" : "s"}`));
  const cuentas = el("div", "cuentas");
  cuentas.appendChild(el("span", "extremo-pep", "H₂N"));
  const abreviatura = {};
  td.codones.forEach((c) => { abreviatura[c.aa] = [c.abrev, c.nombre]; });
  foto.peptido.forEach((aa, i) => {
    if (i > 0) cuentas.appendChild(el("span", "enlace-pep"));
    const [abrev, nombre] = abreviatura[aa];
    const cu = el("span", "cuenta" + (i === foto.peptido.length - 1 && !liberada ? " ultima" : ""), abrev);
    cu.title = nombre;
    cuentas.appendChild(cu);
  });
  cuentas.appendChild(el("span", "extremo-pep", liberada ? "COOH" : "…"));
  cadena.appendChild(cuentas);
  cont.appendChild(cadena);
  return cont;
}

function extraTraduccion(td, foto) {
  const cont = el("div");
  if (td.inicio === -1) return cont;
  const s = foto.sitios;
  const actual = s.A ?? s.P ?? 0;
  const env = el("div", "tabla-envoltura");
  const tabla = el("table", "tabla-codones");
  const cab = el("tr");
  for (const t of ["#", "Codón (5'→3')", "Anticodón (3'→5')", "Aminoácido"]) cab.appendChild(el("th", "", t));
  const thead = el("thead"); thead.appendChild(cab); tabla.appendChild(thead);
  const tbody = el("tbody");
  td.codones.forEach((c, j) => {
    const tr = el("tr", j === (s.P ?? -1) || (c.parada && j === s.A && foto.etapa === "terminacion") ? "actual" : (j > actual ? "pendiente" : ""));
    tr.appendChild(el("td", "", String(j + 1)));
    tr.appendChild(el("td", "mono", c.codon));
    tr.appendChild(el("td", "mono", c.parada ? "(ninguno)" : c.anticodon));
    tr.appendChild(el("td", "", c.parada ? "Parada" : `${c.nombre} (${c.abrev}, ${c.aa})`));
    tbody.appendChild(tr);
  });
  tabla.appendChild(tbody);
  env.appendChild(tabla);
  cont.appendChild(env);
  return cont;
}

function leyendaTraduccion() {
  const lista = el("div", "enzimas-lista");
  const items = [
    ["Ribosoma. ", "Subunidad menor y mayor con tres sitios: A (entra el aminoacil-ARNt), P (ARNt unido a la cadena) y E (sale el ARNt vacío)."],
    ["ARN de transferencia. ", "Su anticodón se empareja con el codón (A-U, C-G) y transporta el aminoácido que corresponde."],
    ["Aminoacil-ARNt sintetasas. ", "Cargan cada ARNt con su aminoácido correcto antes de llegar al ribosoma."],
    ["Factores de liberación. ", "Reconocen los codones de parada UAA, UAG y UGA y liberan la proteína."],
  ];
  for (const [t, d] of items) {
    const div = el("div", "activa");
    div.appendChild(el("strong", "", t));
    div.appendChild(document.createTextNode(d));
    lista.appendChild(div);
  }
  return lista;
}

// ---------------------------------------------------------------------------
// Reproductor
// ---------------------------------------------------------------------------
function fotosEtapa() {
  return estado.datos[estado.etapa].fotos;
}

function mostrar() {
  if (!estado.datos) return;
  const fotos = fotosEtapa();
  const i = Math.min(estado.paso[estado.etapa], fotos.length - 1);
  estado.paso[estado.etapa] = i;
  const foto = fotos[i];

  $("barra").max = String(fotos.length - 1);
  $("barra").value = String(i);
  $("paso-num").textContent = `${i + 1} / ${fotos.length}`;
  $("paso-titulo").textContent = foto.titulo;
  $("paso-texto").textContent = foto.texto;

  const lienzo = $("lienzo");
  const scroll = lienzo.scrollLeft;
  const leyenda = $("leyenda");
  const extra = $("extra");
  lienzo.replaceChildren();
  leyenda.replaceChildren();
  extra.replaceChildren();

  const d = estado.datos;
  if (estado.etapa === "replicacion") {
    lienzo.appendChild(dibujarReplicacion(d.replicacion, foto));
    extra.appendChild(extraReplicacion(d.replicacion, foto));
    leyenda.appendChild(leyendaReplicacion(foto));
  } else if (estado.etapa === "transcripcion") {
    lienzo.appendChild(dibujarTranscripcion(d.transcripcion, foto));
    extra.appendChild(extraTranscripcion(d.transcripcion, foto));
    leyenda.appendChild(leyendaTranscripcion());
  } else {
    lienzo.appendChild(dibujarTraduccion(d.traduccion, foto));
    extra.appendChild(extraTraduccion(d.traduccion, foto));
    leyenda.appendChild(leyendaTraduccion());
  }
  lienzo.scrollLeft = scroll;
  seguirActividad(lienzo);
}

// Desplaza el lienzo para que la zona activa (horquilla, polimerasa, ribosoma) quede a la vista.
function seguirActividad(lienzo) {
  if (!estado.temporizador) return;
  const foco = lienzo.querySelector(".chip.helicasa, .chip.arnpol, .chip.pol3, .codon.sitio-P");
  if (!foco) return;
  const l = lienzo.getBoundingClientRect();
  const f = foco.getBoundingClientRect();
  if (f.left < l.left + 160 || f.right > l.right - 40) {
    lienzo.scrollLeft += f.left - l.left - l.width / 2;
  }
}

function irA(i) {
  const total = fotosEtapa().length;
  estado.paso[estado.etapa] = Math.max(0, Math.min(total - 1, i));
  mostrar();
}

function parar() {
  clearInterval(estado.temporizador);
  estado.temporizador = null;
  const b = document.querySelector('[data-accion="reproducir"]');
  b.textContent = "▶";
  b.setAttribute("aria-label", "Reproducir");
}

function reproducir() {
  if (estado.temporizador) { parar(); return; }
  const total = fotosEtapa().length;
  if (estado.paso[estado.etapa] >= total - 1) estado.paso[estado.etapa] = 0;
  const b = document.querySelector('[data-accion="reproducir"]');
  b.textContent = "⏸";
  b.setAttribute("aria-label", "Pausar");
  estado.temporizador = setInterval(() => {
    const t = fotosEtapa().length;
    if (estado.paso[estado.etapa] >= t - 1) { parar(); return; }
    irA(estado.paso[estado.etapa] + 1);
  }, Number($("velocidad").value));
  mostrar();
}

function cambiarEtapa(etapa) {
  parar();
  estado.etapa = etapa;
  document.querySelectorAll(".etapa-btn").forEach((b) => {
    b.setAttribute("aria-selected", String(b.dataset.etapa === etapa));
  });
  $("lienzo").scrollLeft = 0;
  mostrar();
}

// ---------------------------------------------------------------------------
// Resumen del flujo completo
// ---------------------------------------------------------------------------
function pintarResumen() {
  const d = estado.datos;
  const cont = $("resumen");
  cont.replaceChildren();

  const bloque = (titulo, sub, ...contenido) => {
    const p = el("div", "paso-resumen");
    const r = el("div", "rotulo", titulo);
    r.appendChild(el("small", "", sub));
    const c = el("div");
    contenido.forEach((x) => c.appendChild(x));
    p.append(r, c);
    cont.appendChild(p);
  };
  const linea = (izq, seq, der) => {
    const p = el("div", "seq");
    p.append(document.createTextNode(izq), secuenciaColoreada(seq), document.createTextNode(der));
    return p;
  };

  const rep = d.replicacion, tr = d.transcripcion, td = d.traduccion;
  bloque("ADN de partida", `${d.entrada.longitud} pb, ${d.entrada.gc} % GC`,
    linea("5'-", rep.superior, "-3'"), linea("3'-", rep.inferior, "-5'"));
  bloque("Replicación", "ADN → ADN",
    el("p", "nota", `${rep.fragmentos.length} fragmentos de Okazaki en la cadena rezagada, uno solo continuo en la líder.`),
    el("p", rep.copia_fiel ? "nota ok" : "nota aviso",
      rep.copia_fiel ? "Dos moléculas hijas idénticas a la parental." : "Las copias no coinciden."));
  bloque("Transcripción", `ADN → ARN, molde ${tr.hebra_molde}`, linea("5'-", tr.arnm, "-3'"));

  if (td.inicio === -1) {
    bloque("Traducción", "ARN → proteína",
      el("p", "nota aviso", "El ARNm no tiene ningún AUG, así que no se traduce. Prueba con la otra hebra molde."));
  } else {
    const elementos = [
      el("div", "proteina-final", td.proteina_abrev),
      el("div", "proteina-final", `${td.proteina}  (${td.proteina.length} aa)`),
      el("p", "nota", `Marco de lectura desde el AUG en el nucleótido ${td.inicio + 1} del ARNm.`),
    ];
    if (!td.terminada) elementos.push(el("p", "nota aviso", "No hay codón de parada en el marco: la proteína queda incompleta."));
    bloque("Traducción", "ARN → proteína", ...elementos);
  }
}

// ---------------------------------------------------------------------------
// Entrada y comunicación con el servidor
// ---------------------------------------------------------------------------
const MAXIMO = Number(($("contador").textContent.match(/\/\s*(\d+)/) || [0, 300])[1]);

function actualizarContador() {
  const texto = $("secuencia").value
    .split("\n").filter((l) => !l.trim().startsWith(">")).join("")
    .replace(/[\s\d]/g, "");
  const c = $("contador");
  c.textContent = `${texto.length} / ${MAXIMO} nt`;
  c.classList.toggle("excedido", texto.length > MAXIMO);
}

function mostrarError(msg) {
  const e = $("error");
  e.textContent = msg || "";
  e.hidden = !msg;
}

async function simular() {
  mostrarError("");
  const boton = $("btn-simular");
  boton.disabled = true;
  try {
    const resp = await fetch("/api/simular", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        secuencia: $("secuencia").value,
        tam_fragmento: Number($("tam-fragmento").value),
        tam_cebador: Number($("tam-cebador").value),
        hebra_molde: document.querySelector('input[name="molde"]:checked').value,
      }),
    });
    const datos = await resp.json();
    if (!resp.ok) { mostrarError(datos.error || "No se ha podido simular la secuencia."); return; }
    parar();
    estado.datos = datos;
    estado.paso = { replicacion: 0, transcripcion: 0, traduccion: 0 };
    $("resultado").hidden = false;
    cambiarEtapa("replicacion");
    pintarResumen();
    $("resultado").scrollIntoView({ behavior: matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth" });
  } catch (err) {
    mostrarError("No hay conexión con el servidor. Comprueba que app.py se está ejecutando.");
  } finally {
    boton.disabled = false;
  }
}

async function cargarEjemplos() {
  try {
    const lista = await (await fetch("/api/ejemplos")).json();
    const cont = $("ejemplos");
    lista.forEach((ej, i) => {
      const b = el("button", "ejemplo", ej.nombre);
      b.type = "button";
      b.title = ej.descripcion;
      b.addEventListener("click", () => {
        $("secuencia").value = ej.secuencia;
        actualizarContador();
      });
      cont.appendChild(b);
      if (i === 0 && !$("secuencia").value) { $("secuencia").value = ej.secuencia; actualizarContador(); }
    });
  } catch { /* sin ejemplos: el usuario puede escribir su secuencia */ }
}

function iniciar() {
  $("secuencia").addEventListener("input", actualizarContador);
  $("btn-simular").addEventListener("click", simular);
  $("secuencia").addEventListener("keydown", (e) => {
    if (e.key === "Enter" && (e.ctrlKey || e.metaKey)) simular();
  });
  $("btn-aleatoria").addEventListener("click", async () => {
    try {
      const { secuencia } = await (await fetch("/api/aleatoria")).json();
      $("secuencia").value = secuencia;
      actualizarContador();
    } catch { mostrarError("No se ha podido generar la secuencia."); }
  });

  document.querySelectorAll(".etapa-btn").forEach((b) =>
    b.addEventListener("click", () => cambiarEtapa(b.dataset.etapa)));

  $("reproductor").addEventListener("click", (e) => {
    const accion = e.target.closest("button")?.dataset.accion;
    if (!accion || !estado.datos) return;
    if (accion === "reproducir") { reproducir(); return; }
    parar();
    const i = estado.paso[estado.etapa];
    if (accion === "primero") irA(0);
    if (accion === "anterior") irA(i - 1);
    if (accion === "siguiente") irA(i + 1);
    if (accion === "ultimo") irA(fotosEtapa().length - 1);
  });
  $("barra").addEventListener("input", (e) => { parar(); irA(Number(e.target.value)); });
  $("velocidad").addEventListener("change", () => { if (estado.temporizador) { parar(); reproducir(); } });
  $("lienzo").addEventListener("keydown", (e) => {
    if (!estado.datos) return;
    if (e.key === "ArrowRight") { parar(); irA(estado.paso[estado.etapa] + 1); e.preventDefault(); }
    if (e.key === "ArrowLeft") { parar(); irA(estado.paso[estado.etapa] - 1); e.preventDefault(); }
  });

  cargarEjemplos();
  actualizarContador();
}

iniciar();

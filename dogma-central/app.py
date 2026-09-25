"""Servidor web del simulador del dogma central.

Ejecutar con:  python app.py   y abrir http://127.0.0.1:5000
"""

from flask import Flask, jsonify, render_template, request

from dogma.secuencias import generar_gen_aleatorio
from dogma.simulador import LONGITUD_MAXIMA, simular

app = Flask(__name__)

EJEMPLOS = [
    {
        "nombre": "Inicio de la β-globina humana",
        "descripcion": "Los 11 primeros codones del gen HBB, con un codón de parada añadido.",
        "secuencia": "CTGACGCCACCATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCTAAGCTT",
    },
    {
        "nombre": "β-globina con la mutación falciforme",
        "descripcion": "Mismo fragmento con GAG→GTG en el codón 7: ácido glutámico pasa a valina.",
        "secuencia": "CTGACGCCACCATGGTGCATCTGACTCCTGTGGAGAAGTCTGCCTAAGCTT",
    },
    {
        "nombre": "Gen en la hebra inferior",
        "descripcion": "Con la hebra inferior como molde sale un péptido sin parada; prueba con la superior.",
        "secuencia": "GGATCCTTATTTGCTGGACCAGTAGAACTTCATGGTGGCGAATTC",
    },
]


@app.get("/")
def inicio():
    return render_template("index.html", longitud_maxima=LONGITUD_MAXIMA)


@app.get("/api/ejemplos")
def ejemplos():
    return jsonify(EJEMPLOS)


@app.get("/api/aleatoria")
def aleatoria():
    codones = request.args.get("codones", default=12, type=int)
    return jsonify({"secuencia": generar_gen_aleatorio(codones=codones)})


@app.post("/api/simular")
def api_simular():
    datos = request.get_json(silent=True) or {}
    try:
        resultado = simular(
            datos.get("secuencia", ""),
            tam_fragmento=int(datos.get("tam_fragmento", 8)),
            tam_cebador=int(datos.get("tam_cebador", 3)),
            hebra_molde=datos.get("hebra_molde", "inferior"),
        )
    except (ValueError, TypeError) as e:
        return jsonify({"error": str(e)}), 400
    return jsonify(resultado)


if __name__ == "__main__":
    app.run(debug=True)

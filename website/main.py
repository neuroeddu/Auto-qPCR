from flask import Flask, request, make_response
import io
import readCSV
import pandas as pd


app: Flask = Flask(__name__)


def transform(text_file_contents):
    return text_file_contents.replace("=", ",")


@app.route('/')
def form():
    return ("\n"
            "<html>\n"
            "<head>\n"
            "<meta charset=\"UTF-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">	\n"
            "<title>Untitled Document</title>\n"
            "</head>\n"
            "\n"
            "<body>\n"
            "<link rel=\"stylesheet\"\n"
            "      type=\"text/css\"\n"
            "      href=\"/static/style.css\"/>\n"
            "\n"
            "<header>(Name of the website)</header>\n"                    
            "<table width=\"500\" height=\"250\" border=\"1\">\n"
            "  <tbody>\n"
            "    <tr>\n"
            "      <th width=\"500\" height=\"76\" scope=\"col\">Please upload your data (.csv) </th>\n"
            "    </tr>\n"
            "    <tr>\n"
            "     <td height=\"79\">\n"
            "   <form action=\"/hello\" method=\"post\" enctype=\"multipart/form-data\">\n"    
            "        <input id=\"file-input\" type=\"file\" name= \"data_file\" accept=\".csv\" multiple/><br>\n"
            "	  </td>\n"
            "    </tr>\n"
            "    <tr>\n"
            "      <td height=\"79\">"
            "        Model<br> \n"
            "        <input type=\"radio\" name=\"option\" value=\"hello\" checked>Hello</input><br>\n"
            "        <input type=\"radio\" name=\"option\" value=\"goodbye\" >Goodbye</input><br>\n"
            "        <input type=\"radio\" name=\"option\" value=\"whatsup\" >What's up</input><br>\n"
            "       </td>\n"
            "    </tr>\n"
            "   <tr>\n"
            " <td height=\"79\">"
            "   Name: <br>"
            "   <input type=\"text\" name=\"name\"><br>"
            "       </td>\n"
            "    </tr>\n"
            "   <tr>\n"
            " <td height=\"79\">"
            "   <input type=\"submit\"><br>"
            "</form>"
            "       </td>\n"
            "    </tr>\n"
            "  </tbody>\n"
            "</table>\n"
            "<script>\n"
            "</script>\n"
            "</body>\n"
            "</html>\n")


@app.route('/hello', methods=["POST"])
def transform_view():
    file = request.files['data_file']
    if not file:
        return "No file"

    stream = io.StringIO(file.stream.read().decode("ISO-8859-1"), newline="")
    csv_input = pd.read_csv(stream,
                            skip_blank_lines=True,
                            skipinitialspace=True,
                            engine='python',
                            encoding="utf-8"
                            )
    print(csv_input)
    for i, row in enumerate(csv_input):
        if i <= 1:
            print(row)

    stream.seek(0)
    opt = request.form['option']

    str = request.form['name']
    result = readCSV.csvwriter(stream.read(), opt + " " + str)

    response = make_response(result)
    response.headers["Content-Disposition"] = "attachment; filename=result.csv"
    return response


@app.route('/goodbye', methods=["POST"])
def transform_view2():
    file = request.files['data_file']
    if not file:
        return "No file"

    stream = io.StringIO(file.stream.read().decode("ISO-8859-1"), newline="")
    csv_input = pd.read_csv(stream,
                            skip_blank_lines=True,
                            skipinitialspace=True,
                            engine='python',
                            encoding="utf-8"
                            )
    print(csv_input)
    for i, row in enumerate(csv_input):
        if i <= 1:
            print(row)

    stream.seek(0)
    result = readCSV.csvwriter2(stream.read())

    response = make_response(result)
    response.headers["Content-Disposition"] = "attachment; filename=result.csv"
    return response

from flask import Flask, request, make_response, send_file
import io
import pandas as pd
import AUTOqPCR

app: Flask = Flask(__name__)

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
            "   <form action=\"/download\" method=\"post\" enctype=\"multipart/form-data\">\n"    
            "        <input id=\"file-input\" type=\"file\" name= \"file[]\" accept=\".csv\" multiple/><br>\n"
            "	  </td>\n"
            "    </tr>\n"
            "    <tr>\n"
            "      <td height=\"79\">"
            "        Model<br> \n"
            "        <input type=\"radio\" name=\"option\" value=\"absolute\" checked>Absolute</input><br>\n"
            "        <input type=\"radio\" name=\"option\" value=\"relative\" >Relative</input><br>\n"
            "        <input type=\"radio\" name=\"option\" value=\"stability\" >Stability</input><br>\n"
            "       </td>\n"
            "    </tr>\n"
            "   <tr>\n"
            " <td height=\"79\">"
            "   Control Genes: <br>"
            "   <input type=\"text\" name=\"cgenes\" autocomplete=\"off\"><br>"
            "   Cut-Off: <br>"
            "   <input type=\"text\" name=\"cutoff\" min=\"0\" max=\"1\" step=\"any\" autocomplete=\"off\"><br>"
            "   Max Outliers: <br>"
            "   <input type=\"number\" name=\"max_outliers\" autocomplete=\"off\"><br>"
            "   Control Sample: <br>"
            "   <input type=\"text\" name=\"csample\" autocomplete=\"off\"><br>"
            " </td>\n"
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


@app.route('/download', methods=["POST"])
def transform_view():
    files = request.files.getlist('file[]')
    if len(files) == 0:
        return "No file"

    # Creates empty data frame
    data = pd.DataFrame()

    for item in files:
        stream = io.StringIO(item.stream.read().decode("utf-8") , newline="")
        filedata = pd.read_csv(stream,
                               skip_blank_lines=True ,
                               skipinitialspace=True ,
                               engine='python' ,
                               encoding="utf-8" ,
                               header=46)

        data = data.append(filedata , ignore_index=True , sort=True)
        stream.seek(0)

    model = request.form['option']
    cgenes = request.form['cgenes']
    cutoff = request.form.get('cutoff', type=float)
    max_outliers = request.form.get('max_outliers', type=int)
    csample = request.form['csample']

    result = AUTOqPCR.process_data(data, model, cgenes, cutoff, max_outliers, csample)
    output = result.to_csv()
    response = make_response(output)
    response.headers['Content-Disposition'] = "attachment; filename=output.csv"
    return response


if __name__=='__main__':
    app.debug = True
    app.run(host = '0.0.0.0', port=5000)


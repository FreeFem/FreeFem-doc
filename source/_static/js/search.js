const searchInput = document.getElementById('searchInput')
const searchResults = document.getElementById('searchResults')
let rootPath = ''
const LUNR_LIMIT = 10

function initSearch(path) {
   rootPath = path
}

function search(event) {
   searchClean()

   if (!event.target)
      return

   if (!event.target.value)
      return

   const text = event.target.value

   if (text === '')
      return
   else
      searchLunr(text)
}

function searchClean() {
   searchResults.innerHTML = ''
}

function searchLunr(text) {
   const idx = lunr.Index.load(LUNR_DATA)
   const results = idx.search(text)

   const resultsi = []
   results.forEach(function(result) {
      const ref = result['ref']
      const index = Number(ref)
      const idxi = lunr.Index.load(LUNR_PAGEDATA[index])
      resultsi[index] = idxi.search(text)
   })

   const resultsHTML = parseLunrResults(results, resultsi, text)
   if (resultsHTML.length === 0)
      searchResults.innerHTML = '<p>No results</p>'
   else {
      for (let result of resultsHTML) {
         const div = document.createElement('div')
         div.className = 'search-result'

         const title = document.createElement('a')
         title.className = 'search-result-main-title'
         title.innerHTML = '<i class="fas fa-poll"></i>' + result.mainTitle
         title.href = result.link
         div.appendChild(title)

         for (let i = 0; i < result.titles.length; i++) {
            const subdiv = document.createElement('div')
            subdiv.className = 'search-result-sub'

            const title = document.createElement('a')
            title.className = 'search-result-title'
            title.href = result.links[i]

            const titleDiv = document.createElement('div')

            const titleText = document.createElement('p')
            titleText.className = 'search-result-title-text'
            titleText.innerHTML = result.titles[i]

            const preview = document.createElement('p')
            preview.className = 'search-result-preview'
            preview.innerHTML = result.previews[i]

            title.appendChild(titleDiv)
            titleDiv.appendChild(titleText)
            titleDiv.appendChild(preview)

            subdiv.appendChild(title)
            div.appendChild(subdiv)
         }

         searchResults.appendChild(div)
      }
      searchResults.style.display = 'block'
   }
}

function parseLunrResults(results, resultsi, text) {
   const html = []

   results.forEach(function (result) {
      const ref = result['ref']
      const item = PREVIEW_DATA[ref]
      const mainTitle = item['t']
      const link = rootPath + item['l']

      const titlei = []
      const previewi = []
      const linki = []

      resultsi[ref].forEach(function(subresult) {
         const subref = subresult['ref']
         const subitem = PREVIEW_PAGEDATA[ref][subref]
         const subtitle = subitem['t']
         const sublink = subitem['l']
         const subpreview = subitem['p']

         const lowerCasePreview = subpreview.toLowerCase()
         const lowerCaseText = text.toLowerCase()
         const textIndex = lowerCasePreview.indexOf(lowerCaseText)

         const splittedPreview = subpreview.slice(Math.max(0, textIndex-124), Math.min(textIndex+124, subpreview.length))

         const lowerCasePreview2 = splittedPreview.toLowerCase()
         const textIndex2 = lowerCasePreview2.indexOf(lowerCaseText)
         const finalPreview = splittedPreview.slice(0, textIndex2) + '<b>' + splittedPreview.slice(textIndex2, textIndex2+text.length) + '</b>' + splittedPreview.slice(textIndex2+text.length)

         titlei.push(subtitle)
         previewi.push(finalPreview)
         linki.push(rootPath+sublink)
      })

      html.push({
         mainTitle: mainTitle,
         link: link,
         titles: titlei,
         previews: previewi,
         links: linki
      })
   })



   // for (let i = 0; i < Math.min(LUNR_LIMIT, results.length); i++) {
   //    const id = results[i]['ref']
   //    const item = PREVIEW_DATA[id]
   //    const mainTitle = item['t']
   //    const link = rootPath + item['l']
   //
   //    const titlei = []
   //    const previewi = []
   //    const linki = []
   //    for (let k = 0; k < resultsi.length; k++) {
   //       const results = resultsi[k]
   //       if (results)
   //          for (var j = 0; j < results.length; j++) {
   //             const id = results[j]['ref']
   //             const item = PREVIEW_PAGEDATA[k][id]
   //             const title = item['t']
   //             let preview = item['p']
   //
   //             let lpreview = preview.toLowerCase()
   //             let ltext = text.toLowerCase()
   //             let index = lpreview.indexOf(ltext)
   //             preview = preview.slice(Math.max(0, index-124), Math.min(index+124, preview.length))
   //
   //             lpreview = preview.toLowerCase()
   //             index = lpreview.indexOf(ltext)
   //             preview = preview.slice(0, index) + '<b>' + preview.slice(index, index+text.length) + '</b>' + preview.slice(index+text.length)
   //
   //             const link = rootPath + item['l']
   //             titlei.push(title)
   //             previewi.push(preview)
   //             linki.push(link)
   //          }
   //    }
   //
   //    html.push(
   //       {
   //          mainTitle: mainTitle,
   //          link: link,
   //          titles: titlei,
   //          previews: previewi,
   //          links: linki
   //       }
   //    )
   // }

   return html
}

searchInput.addEventListener('input', function(event) {
   search(event)
})

searchInput.addEventListener('focus', function(event) {
   if (!searchResults.children.length)
     return

   event.stopPropagation()
   searchResults.style.display = 'block'
})

function removeOverlay() {
  document.getElementById("search-overlay").style.display = "none";
}
